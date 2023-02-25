#pragma once

#include <vector>
#include <algorithm>

// https://stackoverflow.com/a/14612943
int sign(int x);

template<typename T> double median(std::vector<T> vec) {
	// This will only be called for small vectors, so just sort.
	std::sort(vec.begin(), vec.end());

	if (vec.size() % 2 == 0) {
		return 0.5 * (vec[vec.size()/2-1] + vec[vec.size()/2]);
	} else {
		return vec[vec.size()/2];
	}
}

// Reads off an unsigned value in most significant byte first format.
unsigned int msb_to_int(const std::vector<unsigned char> & data,
	std::vector<unsigned char>::const_iterator pos, int num_bytes);

unsigned int msb_to_int(const std::vector<unsigned char> & data,
	size_t idx, int num_bytes);

// Calculate a CCITT-CRC16 on the given region.
unsigned short crc16(std::vector<unsigned char>::const_iterator start,
	std::vector<unsigned char>::const_iterator end);
