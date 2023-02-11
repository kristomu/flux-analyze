// Implementation of the Rabin-Karp search algorithm. I might
// add an error-tolerant version later (the principle is easy;
// just generate hashes for each version with errors).

// Worst case complexity is O(nm). Aho-Corasick would be better
// and not restricted to a given string length, but also a much
// greater pain to write.

#include "rabin_karp.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>

// Returns true if there's a needle from pos onwards, false
// otherwise. It assumes the haystack is large enough.
bool rabin_karp::brute_force_search(std::vector<char>::const_iterator pos,
	std::vector<char>::const_iterator end,
	const std::vector<char> & needle) const {

	for (char c: needle) {
		if (pos == end || c != *pos++) {
			return false;
		}
	}

	return true;
}

rabin_karp::rabin_karp(size_t min_needle_length_in) {
	min_needle_length = min_needle_length_in;

	// Calculate the factor required to remove
	// leading (early) characters.
	leading_char_eliminator = 1;

	for (size_t i = 0; i < min_needle_length_in; ++i) {
		leading_char_eliminator = (leading_char_eliminator *
			factor) % modulus;
	}
}

rabin_karp::rabin_karp(const std::vector<char> &
	needle) : rabin_karp(needle.size()) {

	add(needle);
}

void rabin_karp::add(const std::vector<char> & needle) {
	hash_int needle_hash = 0;

	if (needle.size() < min_needle_length) {
		throw std::logic_error("Rabin-Karp: Proposed needle is too short!");
	}

	// If the input needle is longer than the minimum length,
	// then we create a hash that only matches the end; otherwise,
	// we create a hash that matches the whole thing.

	// This works because longer strings are rotated out of the hash
	// anyway, so that at the end of the calculation, only the last
	// min_needle_length bytes are preserved.

	size_t start;
	if (needle.size() > min_needle_length) {
		start = needle.size() - min_needle_length;
	} else {
		start = 0;
	}

	for (size_t i = start; i < needle.size(); ++i) {
		needle_hash = (needle_hash * factor + 
			(unsigned char)needle[i]) % modulus;
	}

	// XXX: If two needles happen to have the same hash,
	// we won't match every needle. Fix later if necessary.

	if (needles_by_hash.find(needle_hash) != needles_by_hash.end()) {
		throw std::invalid_argument("Rabin-Karp: Multiple "
			"strings with same hash not supported!");
	}

	needles_by_hash[needle_hash] = needle;
}

std::vector<size_t> rabin_karp::find_matches(
	const std::vector<char> & haystack) const {

	hash_int search_hash = 0;

	std::vector<size_t> match_indices;

	for (size_t i = 0; i < haystack.size(); ++i) {
		search_hash = (search_hash * factor + 
			(unsigned char)haystack[i]) % modulus;

		// Wraparound might be a problem here... check later
		if (i >= min_needle_length) {
			search_hash -= (leading_char_eliminator * 
				(unsigned char)haystack[i-min_needle_length]) % modulus;

			if (search_hash < 0) {
				search_hash += modulus;
			}
		}

		// If we're far enough that a needle could exist - and we have
		// a hash match, then there might be a match ending in this
		// position.
		if (i + 1 >= min_needle_length &&
			needles_by_hash.find(search_hash) != needles_by_hash.end()) {

			// The +1 is because after the inclusion of the current
			// character above, our "cursor" is now at the character
			// after this one.
			auto needle_ref = needles_by_hash.find(
				search_hash);

			size_t start_pos = i + 1 - needle_ref->second.size();

			if (!brute_force_search(haystack.begin() + start_pos,
				haystack.end(), needle_ref->second)) {
				continue;
			}

			match_indices.push_back(start_pos);
		}
	}

	return match_indices;
}

// Turn a string into a vector, for tests.
std::vector<char> str_to_vec(const std::string & in) {
	return std::vector<char>(in.c_str(), in.c_str() + in.size());
}

template<typename T> bool vec_eq(const T & a, const T & b) {
	if (a.size() != b.size()) { return false; }

	for (size_t i = 0; i < a.size(); ++i) {
		if (a[i] != b[i]) { return false; }
	}

	return true;
}

// I should probably get a unit test framework for this...
void test_rabin_karp() {
	// Test not matching something.
	std::string needle = "isn't here";
	std::string haystack = "hello everybody is this hello world or what? hello";

	rabin_karp rk(str_to_vec(needle));
	if (rk.find_matches(str_to_vec(haystack)).size() > 0) {
		throw std::logic_error("Rabin-Karp test: Found a match "
			"where nowhere was expected");
	}

	// Test matching something.
	needle = "hello";
	rk = rabin_karp(str_to_vec(needle));

	if (!vec_eq(rk.find_matches(str_to_vec(haystack)), {0, 24, 45})) {
		throw std::logic_error("Rabin-Karp test: unexpected output"
			" from three match haystack.");
	}

	// Test needles that are too large to fit in a hash_int.
	needle = "supercalifragilisticexpialidocious";
	haystack = "supercalifragilistic but not "
		"supercalifragilisticexpialidociousexpialidocious";

	rk = rabin_karp(str_to_vec(needle));
	if (!vec_eq(rk.find_matches(str_to_vec(haystack)), {29})) {
		throw std::logic_error("Rabin-Karp test: unexpected output"
			" with very large needle.");
	}

	// Test two needles.
	std::string needletwo = "listic but not";
	rk = rabin_karp(str_to_vec(needletwo));
	rk.add(str_to_vec(needle));

	if (!vec_eq(rk.find_matches(str_to_vec(haystack)), {14, 29})) {
		throw std::logic_error("Rabin-Karp test: unexpected output"
			" with two needles.");
	}

	// Test negative values.
	std::vector<char> haystack_nv = {-1, -1, 0, 1, 1, -1, 0, -1},
		needle_nv = {1, 1, -1, 0, -1};
	rk = rabin_karp(needle_nv);

	if (!vec_eq(rk.find_matches(haystack_nv), {3})) {
		throw std::logic_error("Rabin-Karp test: could not find needle"
			" with negative values.");
	}

}