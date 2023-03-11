#pragma once
#include "search_key.h"
#include <stdexcept>
#include <iostream>

// This is a quick and dirty hybrid of a simple array and a hash
// table with limited double hashing. It's used to store the needles
// (patterns to search for) and optimized for speed, as the main Rabin-
// Karp algorithm will check for a match once every possition in the
// haystack. As we don't need dynamic resizing for our purposes, its
// size is hard-coded.

// The code is thus rather un-C++ like, and closer to ordinary C.
// It should be relatively easy to replace this with a
// std::unordered_map if required for legibility purposes.

class hash_hit {
	public:
		bool empty = true;		// if true, there's nothing here
		hash_int full_hash = 0;
		search_key key;
};

namespace {
	const int MAPSIZE = 65536, DOUBLE_HASH_LIMIT = 10;
}

class search_lookup {
	private:
		hash_hit * table = NULL;

		inline void set_double_hashes(const hash_int & full_hash,
			int & high, int & low) const {

			low = full_hash & (MAPSIZE-1);
			high = (full_hash / MAPSIZE) | 1; // force odd
		}

		inline const hash_hit * get(const int & high,
			const int & low, int double_hash_index) const {

			return &table[(low + double_hash_index * high)
				& (MAPSIZE-1)];
		}

	public:

		search_lookup() {
			table = new hash_hit[MAPSIZE];
		}

		search_lookup & operator=(const search_lookup & other) {
			if (table == NULL) {
				table = new hash_hit[MAPSIZE];
			}
			std::copy(other.table, other.table + MAPSIZE, table);
			return *this;
		}

		search_lookup(const search_lookup & other) {
			*this = other;
		}

		~search_lookup() { delete[] table; }

		void add(hash_int hash_value, search_key key) {
			int high, low;
			set_double_hashes(hash_value, high, low);

			// Find the first empty spot in the hash table.
			int empty_idx;
			for (empty_idx = 0; empty_idx < DOUBLE_HASH_LIMIT &&
				!get(high, low, empty_idx)->empty; ++empty_idx);

			if (empty_idx == DOUBLE_HASH_LIMIT) {
				throw std::logic_error(
					"Hash table is close to full, can't add more");
			}

			int pos = (low + empty_idx * high) & (MAPSIZE-1);
			table[pos].empty = false;
			table[pos].full_hash = hash_value;
			table[pos].key = key;
		}

		const search_key * find(const hash_int & hash_value) const {
			int high, low;
			set_double_hashes(hash_value, high, low);

			for (int i = 0; i < DOUBLE_HASH_LIMIT; ++i) {
				if (get(high, low, i)->empty) {
					return NULL;
				}
				if (get(high, low, i)->full_hash == hash_value) {
					return &get(high, low, i)->key;
				}
			}
			return NULL;
		}
};