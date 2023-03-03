
#include "address_marks.h"
#include "tools.h"

#include <iostream>

// ------------- ID Address Mark ------------

size_t IDAM::byte_length() const { return 10; }

void IDAM::set(std::vector<unsigned char> & raw_bytes) {

	track = raw_bytes[4];
	head = raw_bytes[5];
	sector = raw_bytes[6];

	// Quick and dirty way to deal with out-of-bounds/corrupted values.
	// Also do a lookup to IDAM_datalen to turn this into an actual
	// data chunk length.
	datalen = IDAM_datalen[std::min((int)raw_bytes[7], MAX_DATALEN_IDX)];

	CRC = msb_to_int(raw_bytes, 8, 2);
	CRC_OK = crc16(raw_bytes.begin(), raw_bytes.begin()+8) == CRC;

	truncated = raw_bytes.size() < byte_length();
}

void IDAM::print_info() const {
	std::cout << "ID AM:";
	
	std::cout << " track: " << track << " head: " << head
		<< " sec: " << sector << " data len: " << datalen
		<< " CRC: " << CRC;
	
	if (CRC_OK) {
		std::cout << " OK"; 
	} else {
		std::cout << " bad";
	}
}

// ------------ Data Address Mark (also used for Deleted Data) -----------

size_t DAM::byte_length() const { return 4 + data.size() + 2; }

void DAM::set(std::vector<unsigned char> & raw_bytes, int datalen) {

	deleted = (raw_bytes[3] == 0xF8);
	data = std::vector<char>(raw_bytes.begin() + 4,
		raw_bytes.begin() + 4 + datalen);

	CRC = msb_to_int(raw_bytes, raw_bytes.begin() + 4 + datalen, 2);
	CRC_OK = crc16(raw_bytes.begin(), raw_bytes.begin() + 4 + datalen) == CRC;

	truncated = raw_bytes.size() < byte_length();
}

void DAM::print_info() const {
	if (deleted) {
		std::cout << "Deleted data AM: CRC: " << CRC;	
	} else {
		std::cout << "Data AM: CRC: " << CRC;
	}
	
	if (CRC_OK) {
		std::cout << " OK"; 
	} else {
		std::cout << " bad";
	}
}

// ------------- Index Address Mark ------------

size_t IAM::byte_length() const { return 4; }

void IAM::print_info() const {
	std::cout << "Index AM.";
}

// ------------- Container class ------------

void address_mark::set_address_mark_type(
	const std::vector<unsigned char> & bytes) {

	std::vector<unsigned char> header_bytes(bytes.begin(),
		bytes.begin()+4);

	// The uncommon "WD1772 track language" address marks are from
	// http://dmweb.free.fr/files/Atari-Copy-Protection-V1.4.pdf
	// unless otherwise specified.

	// Remember that these are literal and won't match any corrupted
	// versions. So when we later start reconstructing corrupted
	// address marks, the higher levels of the timeslice must reflect
	// the "pristine" versions, not the corrupted ones.

	mark_type = A_UNKNOWN;
	// Common IAM
	if(header_bytes == std::vector<u_char>({0xC2, 0xC2, 0xC2, 0xFC})) {
		mark_type = A_IAM;
	}
	// Uncommon IAM: http://oldmachinery.blogspot.com/2018/05/qldd.html
	if(header_bytes == std::vector<u_char>({0xC2, 0xC2, 0xC2, 0xFE})) {
		mark_type = A_IAM;
	}

	// Common DDAM
	if(header_bytes == std::vector<u_char>({0xA1, 0xA1, 0xA1, 0xF8})) {
		mark_type = A_DDAM;
	}
	// WD1772 track language: uncommon DDAM
	if(header_bytes == std::vector<u_char>({0xA1, 0xA1, 0xA1, 0xF9})) {
		mark_type = A_DDAM;
	}

	// WD1772 track language: uncommon DAM
	if(header_bytes == std::vector<u_char>({0xA1, 0xA1, 0xA1, 0xFA})) {
		mark_type = A_DAM;
	}
	// Common DAM
	if(header_bytes == std::vector<u_char>({0xA1, 0xA1, 0xA1, 0xFB})) {
		mark_type = A_DAM;
	}

	// Common IDAM
	if(header_bytes == std::vector<u_char>({0xA1, 0xA1, 0xA1, 0xFE})) {
		mark_type = A_IDAM;
	}

	// WD1772 track language: uncommon IDAM
	// Are these real??
	if(header_bytes == std::vector<u_char>({0xA1, 0xA1, 0xA1, 0xFC}) ||
		header_bytes == std::vector<u_char>({0xA1, 0xA1, 0xA1, 0xFD}) ||
		header_bytes == std::vector<u_char>({0xA1, 0xA1, 0xA1, 0xFF})) {
		mark_type = A_IDAM;
	}
}

void address_mark::print_info() const {
	switch(mark_type) {
		default: std::cout << "?? Unknown address mark."; break;
		case A_IAM: return iam.print_info();
		case A_IDAM: return idam.print_info();
		case A_DAM: return dam.print_info();
		case A_DDAM: return ddam.print_info();
	}
}

size_t address_mark::byte_length() const {
	switch(mark_type) {
		default: throw std::runtime_error("Unknown address mark."); break;
		case A_IAM: return iam.byte_length();
		case A_IDAM: return idam.byte_length();
		case A_DAM: return dam.byte_length();
		case A_DDAM: return ddam.byte_length();
	}
}