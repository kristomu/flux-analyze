#include "tools.h"
#include <sstream>

// int to string
std::string itos (int source) {
	std::ostringstream q;
	q << source;
	return (q.str());
}

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