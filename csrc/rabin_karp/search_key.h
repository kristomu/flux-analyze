#pragma once

#include <cstdint>
#include <vector>

// This contains a definition for the hash_int (which is the type
// used to generate "full" hashes), and the search key, which contains
// all the data used to distinguish needles (patterns to find). Each
// needle is associated with a unique ID, so that the code that calls
// Rabin-Karp can find out which needles were found in the haystack.

typedef uint32_t hash_int;

class search_key {
	public:
		std::vector<char> needle;
		int ID;

		search_key(const std::vector<char> & needle_in, int ID_in) {
			needle = needle_in;
			ID = ID_in;
		}

		search_key() {}
};
