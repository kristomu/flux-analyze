--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
----------------------     FLOPPY DISK CONTROLLER    ---------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--                                                                            --
-- Floppy Disk Controller in VHDL.                                            --
-- Copyright (C) <2015>  <Mostafa Abd El-Aziz>                                --
--                                                                            --
-- This program is free software: you can redistribute it and/or modify       --
-- it under the terms of the GNU General Public License as published by       --
-- the Free Software Foundation, either version 3 of the License, or          --
-- (at your option) any later version.                                        --
--                                                                            --
-- This program is distributed in the hope that it will be useful,            --
-- but WITHOUT ANY WARRANTY; without even the implied warranty of             --
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              --
-- GNU General Public License for more details.                               --
--                                                                            --
-- You should have received a copy of the GNU General Public License          --
-- along with this program.  If not, see <http://www.gnu.org/licenses/>.      --
--                                                                            --
-- Contact: iocoder@aol.com                                                   --
--                                                                            --
--------------------------------------------------------------------------------

library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.STD_LOGIC_ARITH.ALL;
use IEEE.STD_LOGIC_UNSIGNED.ALL;

entity fdc is
    Port (
        -- Clock signals
        CLK_50MHz : in  STD_LOGIC;
        CLK_1MHz  : in  STD_LOGIC;

        -- debug
        LED       : out STD_LOGIC_VECTOR (7 downto 0);

        -- CPU Interface
        CS        : in  STD_LOGIC;
        RD        : in  STD_LOGIC;
        WR        : in  STD_LOGIC;
        A         : in  STD_LOGIC_VECTOR (2 downto 0);
        Din       : in  STD_LOGIC_VECTOR (7 downto 0);
        Dout      : out STD_LOGIC_VECTOR (7 downto 0);

        -- Floppy Drive Interface
        DENSEL    : out std_logic := '1'; -- density select
        RESERVED  : in  std_logic;        -- reserved
        INDEX     : in  std_logic;        -- index hole
        MOTEA     : out std_logic := '1'; -- enable motor A
        DRVSB     : out std_logic := '1'; -- select drive B
        DRVSA     : out std_logic := '1'; -- select drive A
        MOTEB     : out std_logic := '1'; -- enable motor B
        DIR       : out std_logic := '1'; -- stepper motor direction
        STEP      : out std_logic := '1'; -- stepper motor pulse
        WDATE     : out std_logic := '1'; -- write data
        WGATE     : out std_logic := '1'; -- write gate
        TRK00     : in  std_logic;        -- track 0
        WPT       : in  std_logic;        -- write protect
        RDATA     : in  std_logic;        -- read data
        SIDE1     : out std_logic := '1'; -- head select
        DSKCHG    : in  std_logic         -- disk change
    );
end fdc;

architecture Behavioral of fdc is

-- read buffer
type   rbuf_t is array (0 to 511) of STD_LOGIC_VECTOR(7 downto 0);
signal rbuf             : rbuf_t                        := (others => x"00");
signal rbuf_ptr         : integer range 0 to 100000     := 0;
signal rbuf_addr        : integer range 0 to 100000     := 0;
signal rbuf_write       : STD_LOGIC                     := '0';
signal rbuf_data        : STD_LOGIC_VECTOR (7 downto 0) := x"00";

-- counter for internal operations:
signal counter          : integer range 0 to 200000000  := 0;

-- phase of execution
signal phase            : integer range 0 to 50         := 0;

-- disk state:
signal motor_on         : std_logic                     := '0';

-- command request sync (if cmd_front != cmd_back, then busy)
signal cmd_front        : integer range 0 to 100000     := 0;
signal cmd_back         : integer range 0 to 100000     := 0;
signal cmd_err          : std_logic                     := '0';

-- read command sync
signal read_front       : integer range 0 to 100000     := 0;
signal read_back        : integer range 0 to 100000     := 0;
signal read_cancel      : STD_LOGIC                     := '0';

-- current track
signal current_track    : integer range 0 to 255        := 0;
signal wanted_track     : integer range 0 to 255        := 0;

-- PLL internal registers:
signal rdata_pulse      : STD_LOGIC                     := '0';
signal last_rdata       : STD_LOGIC                     := '1';
signal length           : integer range 0 to 100000     := 0;
signal delay            : integer range 0 to 100000     := 0;
signal skip             : integer range 0 to 100000     := 0;
signal reset_front      : integer range 0 to 100000     := 0;
signal reset_back       : integer range 0 to 100000     := 0;
signal clock_div        : integer range 0 to 50         := 0;

-- MFM bit buffer
signal mfm_front        : integer range 0 to 100000     := 0;
signal mfm_back         : integer range 0 to 100000     := 0;
signal mfm_type         : integer range 0 to 100000     := 0;
signal mfm_counter      : integer range 0 to 100000     := 0;
signal mfm_clock        : STD_LOGIC                     := '0';
signal mfm_data         : STD_LOGIC                     := '0';

-- read process variables
signal read_active      : STD_LOGIC                     := '0';
signal rd_phase         : integer range 0 to 100000     := 0;
signal rd_subphase      : integer range 0 to 100000     := 0;
signal sync_count       : integer range 0 to 100000     := 0;
signal sync_mark        : STD_LOGIC_VECTOR(47 downto 0) := x"448944894489";
signal sync_indx        : integer range 0 to 50         := 46;
signal rd_indx          : integer range 0 to 50         := 0;
signal is_clock         : STD_LOGIC                     := '1';
signal last_byte        : STD_LOGIC_VECTOR( 7 downto 0) := x"00";
signal tmp_crc          : STD_LOGIC_VECTOR( 7 downto 0) := x"00";
signal crc              : STD_LOGIC_VECTOR(15 downto 0) := x"0000";
signal crc_zeros        : STD_LOGIC                     := '0';

-- FDC registers as seen by CPU:
signal cylinder         : STD_LOGIC_VECTOR (7 downto 0) := x"00";
signal head             : STD_LOGIC_VECTOR (7 downto 0) := x"00";
signal sector           : STD_LOGIC_VECTOR (7 downto 0) := x"00";
signal command          : STD_LOGIC_VECTOR (7 downto 0) := x"00";

begin

--------------------------------------------------------------------------------
--                             FDD Control                                    --
--------------------------------------------------------------------------------

-- translate cylinder and head registers:
wanted_track <= conv_integer(unsigned(cylinder));
SIDE1        <= NOT head(0);

-- translate motor state
MOTEA        <= NOT motor_on;
DRVSA        <= NOT motor_on;

--------------------------------------------------------------------------------
--                            Buffer Access                                   --
--------------------------------------------------------------------------------

process(rbuf_write)
begin
    if (rbuf_write = '0' and rbuf_write'event ) then
        if (read_front /= read_back and read_cancel = '0') then
            rbuf(rbuf_addr) <= rbuf_data;
        end if;
    end if;
end process;

--------------------------------------------------------------------------------
--                           Controller FSM                                   --
--------------------------------------------------------------------------------

process (CLK_50MHz)
begin
    if (CLK_50MHz = '0' and CLK_50MHz'event ) then
        if (phase = 0) then

            -----------------------------
            -- beginning of idle state --
            -----------------------------
            counter  <= 0;
            phase    <= 1;
            STEP     <= '1';
            cmd_back <= cmd_front;

        elsif (phase = 1) then

            -----------------------
            -- idle state itself --
            -----------------------
            if (cmd_front /= cmd_back) then
                -- a command is to be executed.
                counter <= 0;
                if (command(0) = '1') then
                    -- reset
                    phase <= 10;
                elsif (command(1) = '1') then
                    -- read
                    phase <= 20;
                elsif (command(2) = '1') then
                    -- write
                    phase <= 30;
                else
                    -- invalid
                    phase <= 0;
                end if;
            elsif (counter = 100000000) then
                -- it has been idle for more than 2 seconds
                motor_on <= '0';
            else
                counter <= counter + 1;
            end if;

        elsif (phase = 10 or phase = 20 or phase = 30) then

            ---------------
            -- RUN MOTOR --
            ---------------
            if (counter = 0 and motor_on = '1') then
                phase    <= phase + 1;
            elsif (counter = 25000000) then
                counter  <= 0;
                phase    <= phase + 1;
            else
                motor_on <= '1';
                counter  <= counter + 1;
            end if;

        elsif (phase = 11) then

            ----------------------
            -- RESET CONTROLLER --
            ----------------------
            DIR <= '1'; -- from track 1 to track 0
            current_track <= 0;
            if (counter = 0 and TRK00 = '0') then
                phase   <= 0;
            elsif (counter = 25000) then
                STEP    <= '0';
                counter <= counter + 1;
            elsif (counter = 26000) then
                STEP    <= '1';
                counter <= counter + 1;
            elsif (counter = 75000) then
                STEP    <= '0';
                counter <= counter + 1;
            elsif (counter = 76000) then
                STEP    <= '1';
                counter <= counter + 1;
            elsif (counter = 200000) then
                -- seek to next track
                counter <= 0;
            else
                counter <= counter + 1;
            end if;

        elsif (phase = 21 or phase = 31) then

            ----------------
            -- SEEK TRACK --
            ----------------
            if (current_track < wanted_track) then
                DIR <= '0';
            elsif (current_track > wanted_track) then
                DIR <= '1';
            else
                phase <= phase + 1;
            end if;
            if (counter = 25000) then
                STEP    <= '0';
                counter <= counter + 1;
            elsif (counter = 26000) then
                STEP    <= '1';
                counter <= counter + 1;
            elsif (counter = 75000) then
                STEP    <= '0';
                counter <= counter + 1;
            elsif (counter = 76000) then
                STEP    <= '1';
                counter <= counter + 1;
            elsif (counter = 200000) then
                if (current_track < wanted_track) then
                    current_track <= current_track + 1;
                else
                    current_track <= current_track - 1;
                end if;
                counter <= 0;
            else
                counter <= counter + 1;
            end if;

        -----------------
        -- READ SECTOR --
        -----------------
        elsif (phase = 22) then
            phase       <= 23;
            read_cancel <= '0';
            read_front  <= read_front + 1;
        elsif (phase = 23) then
            if (read_front = read_back) then
                -- read operation done
                cmd_err <= '0';
                counter <= 0;
                phase   <= 0;
            elsif (counter = 200000000) then
                -- timeout
                cmd_err     <= '1';
                read_cancel <= '1';
                counter     <= 0;
                phase       <= 0;
            else
                counter <= counter + 1;
            end if;
        end if;


    end if;
end process;

--------------------------------------------------------------------------------
--                      Digital Phase Locked Loop                             --
--------------------------------------------------------------------------------

process(rdata_pulse)

-- length limits
variable one_imp_low  : integer range 0 to 1000 := 0;
variable one_imp_nom  : integer range 0 to 1000 := 0;
variable one_imp_high : integer range 0 to 1000 := 0;

variable two_imp_low  : integer range 0 to 1000 := 0;
variable two_imp_nom  : integer range 0 to 1000 := 0;
variable two_imp_high : integer range 0 to 1000 := 0;

variable tre_imp_low  : integer range 0 to 1000 := 0;
variable tre_imp_nom  : integer range 0 to 1000 := 0;
variable tre_imp_high : integer range 0 to 1000 := 0;

-- delay:
variable old_delay    : integer range 0 to 1000 := 0;
variable new_delay    : integer range 0 to 1000 := 0;

variable length_snap  : integer range 0 to 1000 := 0;

begin

    if (rdata_pulse = '1' and rdata_pulse'event ) then

        if (read_front /= read_back and read_cancel = '0') then

            -- pulse edge, store delay in a variable:
            old_delay    := delay;

            -- calculate limits:
            one_imp_low  := old_delay-old_delay/4;
            one_imp_nom  := old_delay;
            one_imp_high := old_delay+old_delay/4;

            two_imp_low  := old_delay+old_delay/4;
            two_imp_nom  := old_delay+old_delay/2;
            two_imp_high := old_delay+old_delay-old_delay/4;

            tre_imp_low  := old_delay+old_delay-old_delay/4;
            tre_imp_nom  := old_delay+old_delay;
            tre_imp_high := old_delay+old_delay+old_delay/4;

            length_snap  := length;

            -- process MFM bits
            if (skip > 0) then
                skip <= skip - 1;
            elsif (length_snap < delay/4) then
                -- too small
            elsif (length_snap < one_imp_high) then
                -- 1us between two pulses
                new_delay := old_delay - (old_delay - length_snap)/2;
                mfm_type <= 2;
                mfm_front <= mfm_front + 1;
            elsif (length_snap < two_imp_high) then
                -- 1.5us between two pulses
                mfm_type <= 3;
                mfm_front <= mfm_front + 1;
                new_delay := old_delay - (old_delay-(length_snap+length_snap)/4-
                                (length_snap+length_snap)/8 +
                                (length_snap+length_snap)/32)/2;
            elsif (length_snap < delay+delay+delay) then
                -- 2us between two pulses
                mfm_type <= 4;
                mfm_front <= mfm_front + 1;
                new_delay := old_delay - (old_delay - length_snap/2)/2;
            else
                -- too long duration between pulses
                skip <= 5;
            end if;

            -- filter and update delay:
            if (new_delay > old_delay + old_delay/8) then
                delay <= old_delay + old_delay/8;
            elsif (new_delay < old_delay - old_delay/8) then
                delay <= old_delay - old_delay/8;
            else
                delay <= new_delay;
            end if;

            -- initialize length for next pulse
            reset_front <= reset_front + 1;

        else

            -- initialize skip to 1
            skip       <= 1;
            -- initialize delay to 100
            delay      <= 100;

        end if;

    end if;

end process;

process(CLK_50MHz)
begin

    if (CLK_50MHz = '1' and CLK_50MHz'event) then
        if (clock_div = 7) then
            clock_div <= 0;
            rdata_pulse <= rdata;
        else
            if (clock_div = 5) then
                if (reset_front /= reset_back) then
                    length <= 1;
                    reset_back <= reset_front;
                else
                    length <= length + 8;
                end if;
            end if;
            clock_div <= clock_div + 1;
        end if;
    end if;

end process;

--------------------------------------------------------------------------------
--                     Read Command State Machine                             --
--------------------------------------------------------------------------------

process(mfm_clock)
begin
    if (mfm_clock = '0' and mfm_clock'event ) then
        if (read_front /= read_back) then
            if (read_cancel = '1') then
                -- operation cancelled
                read_back   <= read_front; -- finish the deal
                read_active <= '0';
                rbuf_write  <= '0'; -- cancel any write operation
            elsif (read_active = '0') then
                -- operation started
                read_active <= '1';
                rd_phase    <= 1;
            else
                -- currently operating
                if (rd_phase = 1) then
                    ----------------
                    -- INITIALIZE --
                    ----------------
                    sync_indx   <= 46;
                    sync_count  <= 0;
                    rd_subphase <= 0;
                    rd_phase    <= rd_phase + 1;
                elsif (rd_phase = 2) then
                    ----------------------
                    -- LOOKING FOR SYNC --
                    ----------------------
                    if (mfm_data = sync_mark(sync_indx)) then
                        if (sync_indx > 0) then
                            sync_indx  <= sync_indx - 1;
                        else
                            -- sync mark matched
                            rd_phase   <= rd_phase + 1;
                            crc        <= x"E59A"; -- initialize CRC
                            crc_zeros  <= '0';
                            rd_indx    <= 8;   -- will be decremented
                            is_clock   <= '1'; -- next bit is a clock
                        end if;
                    elsif (mfm_data = sync_mark(46)) then
                        sync_indx   <= 45;
                    else
                        sync_indx   <= 46;
                    end if;
                    sync_count <= sync_count + 1;
                elsif (rd_phase = 3) then
                    -------------------------------
                    -- READ SECTOR ADDRESS FIELD --
                    -------------------------------
                    if (is_clock = '1') then
                        -- current bit is a clock
                        is_clock   <= '0';
                        if (rd_indx > 0) then
                            rd_indx <= rd_indx - 1;
                        else
                            rd_indx <= 7;
                            -- *** have just read a full byte ***
                            if (rd_subphase = 0) then
                                -- should be FE
                                if (last_byte = x"FE") then
                                    rd_subphase <= rd_subphase + 1;
                                else
                                    rd_phase  <= 1; -- search again
                                end if;
                            elsif (rd_subphase = 1) then
                                -- track
                                if (last_byte = cylinder) then
                                    rd_subphase <= rd_subphase + 1;
                                else
                                    rd_phase  <= 1; -- search again
                                end if;
                            elsif (rd_subphase = 2) then
                                -- head
                                if (last_byte = head) then
                                    rd_subphase <= rd_subphase + 1;
                                else
                                    rd_phase  <= 1; -- search again
                                end if;
                            elsif (rd_subphase = 3) then
                                -- sector
                                if (last_byte = sector) then
                                    rd_subphase <= rd_subphase + 1;
                                else
                                    rd_phase  <= 1; -- search again
                                end if;
                            elsif (rd_subphase = 4) then
                                -- sector length
                                if (last_byte = x"02") then
                                    rd_subphase <= rd_subphase + 1;
                                    crc_zeros   <= '1';
                                else
                                    rd_phase  <= 1; -- search again
                                end if;
                            elsif (rd_subphase = 5) then
                                -- CRC 1
                                tmp_crc <= last_byte;
                                rd_subphase <= rd_subphase + 1;
                            elsif (rd_subphase = 6) then
                                -- CRC 2
                                if ((tmp_crc & last_byte) = crc) then
                                    rd_subphase <= rd_subphase + 1;
                                else
                                    rd_phase  <= 1; -- search again
                                end if;
                            elsif (rd_subphase = 7) then
                                -- seek to next A1 address mark
                                sync_indx   <= 46;
                                sync_count  <= 0;
                                rd_phase    <= 2;
                                rd_subphase <= rd_subphase + 1;
                            elsif (rd_subphase = 8) then
                                if (sync_count > 700) then
                                    rd_phase  <= 1; -- search again
                                elsif (last_byte=x"FB" or last_byte=x"F8") then
                                    rd_subphase <= rd_subphase + 1;
                                else
                                    rd_phase  <= 1; -- search again
                                end if;
                            elsif (rd_subphase > 8 and rd_subphase < 521) then
                                rbuf_addr   <= rd_subphase-9;
                                rbuf_data   <= last_byte;
                                rbuf_write  <= '1';
                                if (rd_subphase = 520) then
                                    crc_zeros <= '1';
                                end if;
                                rd_subphase <= rd_subphase + 1;
                            elsif (rd_subphase = 521) then
                                -- CRC 1
                                tmp_crc     <= last_byte;
                                rd_subphase <= rd_subphase + 1;
                            elsif (rd_subphase = 522) then
                                -- CRC 2
                                if ((tmp_crc & last_byte) = crc) then
                                    rd_subphase <= rd_subphase + 1;
                                else
                                    rd_phase  <= 1; -- search again
                                end if;
                            elsif (rd_subphase = 523) then
                                -- finally done
                                read_back   <= read_front;
                                read_active <= '0';
                            end if;
                        end if;
                    else
                        -- read data bit
                        last_byte(rd_indx) <= mfm_data;
                        -- calculate CRC
                        if (crc_zeros = '0') then
                            if (crc(15) = '1') then
                                crc <= (crc(14 downto 0)&mfm_data) XOR x"1021";
                            else
                                crc <= (crc(14 downto 0)&mfm_data);
                            end if;
                        else
                            if (crc(15) = '1') then
                                crc <= (crc(14 downto 0)&"0") XOR x"1021";
                            else
                                crc <= (crc(14 downto 0)&"0");
                            end if;
                        end if;
                        -- next bit is a clock
                        is_clock   <= '1';
                        -- stop any writing to rbuf
                        rbuf_write <= '0';
                    end if;
                end if;
            end if;
        end if;
    end if;
end process;

process(CLK_50MHz)
begin
    if (CLK_50MHz = '1' and CLK_50MHz'event ) then
        if (mfm_back /= mfm_front) then
            if (mfm_counter = 0) then
                mfm_data  <= '1';
                mfm_clock <= '1';
            elsif (mfm_counter = 1) then
                mfm_clock <= '0';
            elsif (mfm_counter = 2) then
                mfm_data  <= '0';
                mfm_clock <= '1';
            elsif (mfm_counter = 3) then
                mfm_clock <= '0';
            elsif (mfm_counter = 4) then
                mfm_clock <= '1';
            elsif (mfm_counter = 5) then
                mfm_clock <= '0';
            elsif (mfm_counter = 6) then
                mfm_clock <= '1';
            elsif (mfm_counter = 7) then
                mfm_clock <= '0';
            end if;
            if (mfm_counter = mfm_type*2 - 1) then
                mfm_back    <= mfm_front;
                mfm_counter <= 0;
            else
                mfm_counter <= mfm_counter + 1;
            end if;
        else
            mfm_counter <= 0;
            mfm_clock   <= '0';
        end if;
    end if;
end process;

--------------------------------------------------------------------------------
--                           CPU Interface                                    --
--------------------------------------------------------------------------------

process (CLK_50MHz)
begin
    if (CLK_50MHz = '0' and CLK_50MHz'event) then
        if (CS = '1') then
            if (A = "000" and WR = '1') then
                -- cylinder register
                cylinder <= Din;
            elsif (A = "001" and WR = '1') then
                -- head register
                head     <= Din;
            elsif (A = "010" and WR = '1') then
                -- sector register
                sector   <= Din;
            elsif (A = "011" and WR = '1') then
                -- command register
                command   <= Din;
                cmd_front <= cmd_front + 1;
                rbuf_ptr  <= 0;
            elsif (A = "100" and RD = '1') then
                -- status register
                if (cmd_front = cmd_back) then
                    Dout(0) <= '0'; -- not busy
                else
                    Dout(0) <= '1'; -- busy
                end if;
                Dout(1)          <= cmd_err;
                Dout(7 downto 2) <= "000000";
            elsif (A = "101" and RD = '1') then
                -- data out
                Dout <= rbuf(rbuf_ptr);
            elsif (A = "110" and WR = '1') then
                -- data in
            elsif (A = "111" and WR = '1') then
                -- update pointer
                if (rbuf_ptr = 511) then
                    rbuf_ptr <= 0;
                else
                    rbuf_ptr <= rbuf_ptr + 1;
                end if;
            end if;
        end if;
    end if;
end process;

end Behavioral;