// Standalone FluxEngine IBM decoder with some tricks
// from the Python version... hopefully faster and
// more accurate. Perhaps multiple flux file support ???

#include <iostream>
#include <sqlite3.h>

#include <string>
#include <stdexcept>
#include <vector>
#include <cmath>

#include <zlib.h>

// Taken from a zlib example. If I really want to optimize stuff,
// there exist faster implementations of inflate than zlib's.
std::string decompress_zlib(const std::string & input) {

	// Initialize zstream
	z_stream zstream;

	zstream.opaque = Z_NULL;
	zstream.zalloc = Z_NULL;
	zstream.zfree = Z_NULL;

	zstream.avail_in = 0;
	zstream.next_in = Z_NULL;

	zstream.avail_out = 0; // pump priming

	int ret = inflateInit(&zstream);
	if (ret != Z_OK) {
		throw std::runtime_error("Could not initialize zlib stream! Error: " +
			std::string(zError(ret)));
	}

	// Set the input parameters.
	zstream.next_in = (Bytef *)input.c_str();
	zstream.avail_in = input.size();

	std::string decompressed;

	// Decompress a buffer's worth at a time.
	const int BUFLEN = 32768;

	char decompressed_buf[BUFLEN];
	while (zstream.avail_out == 0) {
		zstream.avail_out = BUFLEN;
		zstream.next_out = (Bytef *) decompressed_buf;

		ret = inflate(&zstream, Z_NO_FLUSH);

		if (ret != Z_OK && ret != Z_STREAM_END) {
			inflateEnd(&zstream);
			throw std::runtime_error("Could not decompress data! Error: " +
				std::string(zError(ret)));
		}

		decompressed += std::string(decompressed_buf,
			decompressed_buf + BUFLEN - zstream.avail_out);
	}

	// All done, do some cleanup, and then return.
	inflateEnd(&zstream);
	return decompressed;
}

enum comp_type { CT_UNCOMPRESSED = 0, CT_ZLIB = 1 };

class flux_data {
	public:
		int track, side;
		std::vector<int> flux_delays;

		flux_data(int track_in, int side_in,
			comp_type compression_type, std::string fluxengine_data) {

			track = track_in;
			side = side_in;

			std::string uncompressed_data;

			// Decompress if necessary.
			switch(compression_type) {
				case CT_UNCOMPRESSED:
					uncompressed_data = fluxengine_data;
					break;
				case CT_ZLIB:
					uncompressed_data = decompress_zlib(fluxengine_data);
					break;
				default:
					// TODO: itos
					throw std::runtime_error(
						"Unknown FluxEngine compression type!");
			}

			// As far as I understand it, the FluxEngine native format
			// is simply: take each byte & 0x3F. If it is 0x3F, then that
			// just indicates that the flux delay is larger than 63, and
			// so we should accumulate them. Otherwise the flux delay
			// is just the byte & 0x3F.
			
			flux_delays.reserve(uncompressed_data.size());

			int next_delay = 0;
			for (unsigned char data_byte: uncompressed_data) {
				next_delay += data_byte & 0x3F;
				if ((data_byte & 0x3F) != 0x3F) {
					flux_delays.push_back(next_delay);
					next_delay = 0;
				}
			}

			if (next_delay != 0) {
				flux_delays.push_back(next_delay);
			}
		}
};

// https://stackoverflow.com/a/14612943
int sign(int x) {
	return (x > 0) - (x < 0);
}

// Higher order formatting for IBM floppies, from
// https://www-user.tu-chemnitz.de/~heha/basteln/PC/usbfloppy/floppy.chm/

// The start of a chunk is marked by a synchronization header. There are
// four types:

// Abbrev. Full name                 Byte seq.  What it is
// IAM	   Index Address Mark        C2C2C2 FC, start of a track
// IDAM	   ID Address Mark           A1A1A1 FE, start of a sector
// DAM	   Data Address Mark         A1A1A1 FB, start of sector data
//		   Deleted Data Address Mark A1A1A1 F8,	start of deleted sector data

// Each address mark is followed by a different header structure,
// though the IAM has no header at all. Note that the byte sequence
// contains a deliberate error to signal out-of-band that it's an address
// mark, not ordinary data.

// Because we want to be robust, we're going to only search for the short
// MFM sequences; we don't make any assumptions about the gap data that
// precedes them.

// The following are given as bit sequences, i.e. *before* MFM decoding.
// Hence we don't have to deal with detecting out-of-band errors.

// // //

// Turn a sequence (e.g. 1 0 1 0 0 1) into a run length list (1 2). Run
// length lists are used for the ordinal search approach.

// I really should describe how the system is designed...

std::vector<int> sequence_to_run_lengths(const std::vector<char> & sequence) {
	// The sequence must start and end in an 1, because otherwise
	// we don't know if, say "0 1" isn't really an excerpt from
	// "0 0 1"; ditto, we don't know if "1 0" isn't really an excerpt
	// from "1 0 0 ". So first establish this.

	//if (sequence.size() == 0) { return {}; }

	size_t begin, end;

	for (begin = 0; begin < sequence.size() &&
		sequence[begin] != 1; ++begin);

	for (end = sequence.size()-1; end+1 > 0 &&
		sequence[end] != 1; --end);

	std::vector<int> run_lengths;

	int cur_run_length = 0;

	for (size_t i = begin+1; i <= end; ++i) {
		if (sequence[i] == 0) {
			++cur_run_length;
		} else {
			run_lengths.push_back(cur_run_length);
			cur_run_length = 0;
		}
	}

	return run_lengths;
}

class MFM_magic_bytes {
	public:
		std::vector<char> short_A1 = {
			0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1};
		std::vector<char> short_C2 = {
			0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0};
		std::vector<char> A1_sequence, C2_sequence;
		std::vector<int> A1_run_lengths, C2_run_lengths;

		MFM_magic_bytes();
};

MFM_magic_bytes::MFM_magic_bytes() {
	// The A1 sequence is the short A1 repeated three times,
	// and ditto for the C2 sequence.

	A1_sequence.reserve(3 * short_A1.size());
	C2_sequence.reserve(3 * short_C2.size());
	for (int i = 0; i < 3; ++i) {
		std::copy(short_A1.begin(), short_A1.end(),
			std::back_inserter(A1_sequence));
		std::copy(short_C2.begin(), short_C2.end(),
			std::back_inserter(C2_sequence));
	}

	A1_run_lengths = sequence_to_run_lengths(A1_sequence);
	C2_run_lengths = sequence_to_run_lengths(C2_sequence);
}

void ordinal_search(const flux_data & flux) {

	// First generate a delta array, which consists of the sign value
	// of adjacent values of the flux delay. If a delay is increasing,
	// then the delta is +1; if it's decreasing, it's -1; if it's the
	// same, it's 0.
	
	std::vector<int> flux_delta(flux.flux_delays.size()-1, 0);

	for (size_t i = 0; i < flux_delta.size(); ++i) {
		flux_delta[i] = sign(flux.flux_delays[i+1] - flux.flux_delays[i]);
		std::cout << flux_delta[i] << " ";
	}

	std::cout << std::endl;
}

// PROCESSING

// There are two stage from here. The first is just the train of reversals/
// no flux reversals, e.g. if the clock is 20 and we have 40, then we first
// have no reversal (0) then a reversal (1). In a way this is the "integral"
// of the flux deltas, quantized wrt. the clock.
// Unlike the Python program, I won't be using three bands, I'll just be using
// one clock; this is the easiest way to do things. I can test whether there's
// a benefit to going to the more complex three-band version later.
// Note also that the length of the first clock is undefined. E.g. if the first
// byte is 40, then that means that, from the point where the read head started
// reading to the first flux transition, two clocks happened. But we have no
// way of knowing where the previous flux transition happened.

// The second is the MFM coding itself. Consider the bit pattern that the
// first stage produces; then
// https://info-coach.fr/atari/hardware/FD-Hard.php#Modified_Frequency_Modulation
// has a good table of how to turn this into MFM.
// If there's an NN, i.e. 00, and the last bit was 0, then that's an error.
// It's used deliberately in sync sequences but otherwise indicates that
// something went wrong. The same goes for RR, which I've (arbitrarily?)
// designated a one; it should never happen. There may be better ways to
// handle RRs, e.g. assuming everything shorter than the clock is RNR.

// Error_out is an output parameter where the error (departure from perfect
// fit) is stored.
//
// We should handle very long reversals, e.g. RNNNNNNR, the expected way
// even though IBM floppies can't have more than three non-reversals between
// each reversal; this because I may want to overlap different reads of the
// same track later to combine them to find otherwise missing reversals. Note
// that from an ideal perspective, a flux reversal is infinitesimally short
// because the square wave it's differentiated from would be perfectly sharp.

// TODO: When searching for patterns (i.e. not correcting errors), it may be
// quicker to make a run-length encoding of this (e.g. 10010001 = 2 3, and then
// search for run-length encoded stuff instead, since both the haystack and
// needle will be shorter. Do that perhaps?)

std::vector<char> get_sequence_bits(double clock,
		const std::vector<int> delays, double & error_out) {

	std::vector<char> sequence_bits;

	error_out = 0;

	// Skip the first delay byte (see above). We know there's a
	// flux transition after the first delay, so our stream begins
	// with a one.
	sequence_bits.reserve(delays.size() * 4 / 3); // expected number of bits per reversal
	sequence_bits.push_back(1);

	for (size_t i = 1; i < delays.size(); ++i) {
		double frac_clock = delays[i]/clock;
		int period = std::round(frac_clock);
		double error_term = delays[i] - period * clock;

		//std::cout << "PARSE TEST " << (int)delays[i] << " " << period << std::endl;

		// PLL test.
		//clock = 0.9 * clock + 0.1 * delays[i]/(double)period;
		//std::cout << clock << std::endl;

		// Update Euclidean error.
		error_out += error_term * error_term;

		for (int j = 0; j < period; ++j) {
			sequence_bits.push_back(0);
		}
		sequence_bits.push_back(1);
	}

	error_out = std::sqrt(error_out / delays.size());

	return sequence_bits;
}

// Read the flux data from the given FluxEngine filename. Returns a
// vector of flux data. If verbose is set to true, output stats as
// the data is read.

std::vector<flux_data> get_flux_data(std::string flux_filename,
	bool verbose) {

	sqlite3 * database;

	// Open FluxEngine sqlite database.
	if (sqlite3_open(flux_filename.c_str(), &database) != SQLITE_OK) {
		sqlite3_close(database);
		throw std::runtime_error("Could not open flux file " + flux_filename);
	}

	// Now get every track and side. (NOTE: This differs from the MFM blocks,
	// which are categorized by sector, IIRC.)
	sqlite3_stmt * statement;
	std::string flux_query = "SELECT track, side, data, compression FROM zdata";
	if (sqlite3_prepare_v2(database, flux_query.c_str(),
		-1, &statement, NULL) != SQLITE_OK) {

		std::string err = "Could not prepare statement! sqlite error: " +
			std::string(sqlite3_errmsg(database));
		sqlite3_close(database);

		throw std::runtime_error(err);
	}

	std::vector<flux_data> flux_records;

	bool done = false;
	while (!done) {
		switch(sqlite3_step(statement)) {
			case SQLITE_DONE:
				done = true;
				break;
			case SQLITE_ROW: {
				int track = sqlite3_column_int(statement, 0);
				int side = sqlite3_column_int(statement, 1);
				comp_type compression = (comp_type)
					sqlite3_column_int(statement, 3);

				int data_size = sqlite3_column_bytes(statement, 2);
				const char * data = (const char *)sqlite3_column_blob(
					statement, 2);
				std::string data_str(data, data + data_size);

				flux_data this_row(track, side, compression, data_str);

				if (verbose) {
					std::cout << "T: " << track << ", S: "
						<< side << ", compression: " << compression
						<< ", num bytes: " << data_size << ", total bytes: " <<
						this_row.flux_delays.size() << std::endl;
				}

				// TODO: Stuff it into a proper data structure and handle it elsewhere.

				flux_records.push_back(this_row);
				break;
			}
			default: {
				std::cout << sqlite3_step(statement) << std::endl;
				std::string err = "Could not retrieve row! sqlite error: " +
					std::string(sqlite3_errmsg(database));
				sqlite3_close(database);

				throw std::runtime_error(err);
			}
		}
	}

	// All done, get outta here.
	sqlite3_finalize(statement);
	sqlite3_close(database);

	return flux_records;
}

/* And then the plan is something like:
 * 	- Keep a list of MFM intro markers
 *	- Generate ordinal delta patterns
 *	- Use these as needles and the flux delta as haystack
 *	- For each match, do a narrowing down to determine the clock 
 *	- Try to decode the match. If OK, all's well. Otherwise... */

/* Then later add fancy stuff like warping (PLL), possible lining up
 * multiple incomplete reads - e.g. if #1 is 2 2 8 and #2 is 2 5 4, then
 * #1 = ..X..X.......X
 * #2 = ..X.....X....X
 * which could suggest
 *    = ..X..X..X....X
 * as a true composition. */

int vmain() {
	std::string flux_to_read = "VIGOR.flux"; // fix later

	sqlite3 * database;

	// Open FluxEngine sqlite database.
	if (sqlite3_open(flux_to_read.c_str(), &database) != SQLITE_OK) {
		sqlite3_close(database);
		throw std::runtime_error("Could not open flux file " + flux_to_read);
	}

	// Now get every track and side. (NOTE: This differs from the MFM blocks,
	// which are categorized by sector, IIRC.)
	sqlite3_stmt * statement;
	std::string flux_query = "SELECT track, side, data, compression FROM zdata";
	if (sqlite3_prepare_v2(database, flux_query.c_str(),
		-1, &statement, NULL) != SQLITE_OK) {

		std::string err = "Could not prepare statement! sqlite error: " +
			std::string(sqlite3_errmsg(database));
		sqlite3_close(database);

		throw std::runtime_error(err);
	}

	std::vector<flux_data> flux_records;

	bool done = false;
	while (!done) {
		switch(sqlite3_step(statement)) {
			case SQLITE_DONE:
				done = true;
				break;
			case SQLITE_ROW: {
				int track = sqlite3_column_int(statement, 0);
				int side = sqlite3_column_int(statement, 1);
				comp_type compression = (comp_type)
					sqlite3_column_int(statement, 3);

				int data_size = sqlite3_column_bytes(statement, 2);
				const char * data = (const char *)sqlite3_column_blob(
					statement, 2);
				std::string data_str(data, data + data_size);

				std::cout << "T: " << track << ", S: " << side << ", compression: "
					<< compression << ", num bytes: " << data_size;

				// TODO: Stuff it into a proper data structure and handle it elsewhere.
				flux_data foo(track, side, compression, data_str);
				std::cout << std::endl << "\t" << foo.flux_delays.size() << std::endl;
				flux_records.push_back(foo);
				break;
			}
			default: {
				std::cout << sqlite3_step(statement) << std::endl;
				std::string err = "Could not retrieve row! sqlite error: " +
					std::string(sqlite3_errmsg(database));
				sqlite3_close(database);

				throw std::runtime_error(err);
			}
		}
	}

	// All done, get outta here.
	sqlite3_finalize(statement);
	sqlite3_close(database);

	// Stuff not related to the database goes here.
	//
	for (const flux_data & f: flux_records) {
		double error;

		std::vector<char> sequence = get_sequence_bits(23,
				f.flux_delays, error);

		std::cout << f.track << ", " << f.side << " error = " << error << "\n";
	}
	return 0;
}

// Quick and dirty (brute force) string search. Fix later.
std::vector<int> search(const std::vector<char> & haystack,
	const std::vector<char> & needle) {

	std::vector<int> hits;

	for (size_t i = 0; i < haystack.size(); ++i) {
		std::cerr << i / (double)haystack.size() << "   \r" << std::flush;
		bool match = true;
		for (size_t j = 0; j < needle.size() && match; ++j) {
			match &= haystack[i+j] == needle[j];
		}
		if (match) {
			hits.push_back(i);
		}
	}

	return hits;
}

int main() {

	std::vector<flux_data> flux_records = get_flux_data(
		"../tracks/low_level_format_with_noise.flux", true);
	// Stuff not related to the database goes here.

	MFM_magic_bytes magic;

	for (const flux_data & f: flux_records) {
		double error;

		std::vector<char> sequence = get_sequence_bits(17,
				f.flux_delays, error);

		std::cout << f.track << ", " << f.side << " error = " << error << "\n";

		// Do a simple search for A1 patterns.
		// Currently inordinately slow due to O(mn) complexity. I am aware :-)
		std::vector<int> a1_positions = search(sequence, magic.A1_sequence);
		std::cout << "Sequences found: " << a1_positions.size() << std::endl;
		std::cout << std::endl;

	}
	return 0;
}