// Implementation of the Rabin-Karp search algorithm. I might
// add an error-tolerant version later (the principle is easy;
// just generate hashes for each version with errors).

// Worst case complexity is O(nm). Aho-Corasick would be better
// and not restricted to a given string length, but also a much
// greater pain to write.

#pragma once
#include <vector>
#include <iostream>
#include <cstdint>
#include <stdexcept>

#include "hashtable.h"

typedef uint32_t hash_int;

// Each needle is associated with a different ID so that
// it's possible to tell just what string was matched.

class search_result {
	public:
		size_t idx;
		int ID;

		search_result(size_t idx_in, int ID_in) {
			idx = idx_in;
			ID = ID_in;
		}

		bool operator==(const search_result & x) const {
			return idx == x.idx && ID == x.ID;
		}

		bool operator!=(const search_result & x) const {
			return !(*this == x);
		}

		search_result() {}
};

class rabin_karp {
	private:
		search_lookup needles_by_hash;

		// This is the magic constant of the Bernstein hash. See
		// https://azrael.digipen.edu/~mmead/www/Courses/CS280/HashFunctions-1.html
		// for discussion.
		static const hash_int factor = 33;

		size_t min_needle_length;
		hash_int leading_char_eliminator;

		bool brute_force_search(std::vector<char>::const_iterator pos,
			std::vector<char>::const_iterator end,
			const std::vector<char> & needle) const;

	public:
		// Each needle is associated with an ID, and the find_matches
		// function returns search results that contain the ID and
		// index (location) of the match.

		rabin_karp(size_t min_needle_length_in);
		rabin_karp(const std::vector<char> & needle, int needle_ID);
		void add(const std::vector<char> & needle, int needle_ID);
		std::vector<search_result> find_matches(
			const std::vector<char> & haystack,
			size_t max_matches) const;
		std::vector<search_result> find_matches(
			const std::vector<char> & haystack) const;
		std::vector<size_t> find_matching_indices(
			const std::vector<char> & haystack) const;
};

// Unit test. I should probably get a proper test framework for this...
void test_rabin_karp();