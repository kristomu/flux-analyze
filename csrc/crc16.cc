#include "crc16.h"
#include <stdexcept>
#include <iostream>

std::vector<unsigned short> crc16_table = create_crc16_table();

unsigned short crc16_bitwise(
	std::vector<unsigned char>::const_iterator start,
	std::vector<unsigned char>::const_iterator end,
	bool creating_table) {

	short i, crc = -1;

	if (creating_table) {
		// If we're creating a table, the initial CRC value
		// must be 0 since we can't chain the one-byte CRCs
		// together otherwise. We'll set the initial value
		// inside the table-based CRC function instead.
		crc = 0;
	}

	for (auto pos = start; pos != end; ++pos) {
		crc ^= (int)*pos << 8;
		i = 8;
		do {
			if (crc & 0x8000) {
				crc = crc << 1 ^ 0x1021;
			} else {
				crc = crc << 1;
			}
		} while (--i);
	}

	return crc;
}

unsigned short crc16_bitwise(
	std::vector<unsigned char>::const_iterator start,
	std::vector<unsigned char>::const_iterator end) {

	return crc16_bitwise(start, end, false);
}

std::vector<unsigned short> create_crc16_table() {
	// Cache the CRCs of each byte value so we can just look them up
	// instead of processing per bit. The same idea can be extended
	// to 16 bits if required, but I don't think that's necessary.

	std::vector<unsigned short> table;
	std::vector<unsigned char> char_to_CRC(1);

	for (int i = 0; i < 256; ++i) {
		char_to_CRC[0] = i;
		table.push_back(crc16_bitwise(char_to_CRC.begin(),
			char_to_CRC.end(), true));
	}

	return table;
}

unsigned short crc16(
	std::vector<unsigned char>::const_iterator start,
	std::vector<unsigned char>::const_iterator end) {

	if (crc16_table.empty()) {
		crc16_table = create_crc16_table();
	}

	unsigned short crc = -1;

	for (auto pos = start; pos != end; ++pos) {
		crc = (crc << 8) ^ crc16_table[((crc >> 8) ^ *pos)];
	}

	return crc;
}

void test_crc16() {
	std::vector<unsigned char> test = {'t', 'e', 's', 't'};

	unsigned short table_crc = crc16(test.begin(), test.end());
	unsigned short bitwise_crc = crc16_bitwise(test.begin(),
		test.end());

	if (table_crc != bitwise_crc) {
		throw std::logic_error("CRC16 test: Discrepancy between"
			" table-based and bitwise calculation!");
	}
}