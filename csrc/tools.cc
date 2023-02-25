#include "tools.h"

int sign(int x) {
	return (x > 0) - (x < 0);
}

// Reads off an unsigned value in most significant byte first format.
unsigned int msb_to_int(const std::vector<unsigned char> & data,
	std::vector<unsigned char>::const_iterator pos, int num_bytes) {

	unsigned int out = 0;

	for (int i = 0; i < num_bytes; ++i) {
		if (pos == data.end()) { return out; } // early stopping
		out <<= 8;
		out += (int)*pos++;
	}

	return out;
}

unsigned int msb_to_int(const std::vector<unsigned char> & data,
	size_t idx, int num_bytes) {
	return msb_to_int(data, data.begin()+idx, num_bytes);
}

// https://stackoverflow.com/a/44947877/6183577
// But as this is CCITT-CRC16, the initial CRC value is -1 (0xFFFF), not 0.
unsigned short crc16(std::vector<unsigned char>::const_iterator start,
	std::vector<unsigned char>::const_iterator end) {

	short i, crc = -1;

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