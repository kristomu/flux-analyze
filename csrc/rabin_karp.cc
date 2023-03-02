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
	needle, int needle_ID) : rabin_karp(needle.size()) {

	add(needle, needle_ID);
}

void rabin_karp::add(const std::vector<char> & needle, int needle_ID) {
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

	needles_by_hash[needle_hash] = search_key(needle, needle_ID);
}

std::vector<search_result> rabin_karp::find_matches(
	const std::vector<char> & haystack) const {

	hash_int search_hash = 0;

	std::vector<search_result> matches;

	// The boundary conditions for this loop are tricky. When we
	// enter the loop with some value of i, that means that every
	// byte before i has been hashed. So a needle whose last byte
	// is at position i-1 will now be detected - this needle has
	// length j and starts at i-j. (E.g. length 8, detected at i=8,
	// starts at 0.) Hence we need one extra round through the loop
	// to detect any needles at exactly the end of the haystack.

	auto needle_ref = needles_by_hash.find(0);

	for (size_t i = 0; i <= haystack.size(); ++i) {
		// If we're far enough that a needle could exist - and we have
		// a hash match, then there might be a match ending in this
		// position.
		needle_ref = needles_by_hash.find(search_hash);

		// The possible match needs to be longer than the minimum length,
		// it needs to have a hash that we know,
		// and it needs to be longer than the needle for that hash.

		if (i >= min_needle_length &&
			needle_ref != needles_by_hash.end() &&
			needle_ref->second.needle.size() <= i) {

			size_t start_pos = i - needle_ref->second.needle.size();

			// Eliminate false positives by verifying with brute force.
			if (brute_force_search(haystack.begin() + start_pos,
				haystack.end(), needle_ref->second.needle)) {

				matches.push_back(search_result(start_pos,
					needle_ref->second.ID));
			}
		}

		// If we're at the end of the haystack, don't add any more
		// bytes.
		if (i == haystack.size()) { continue; }

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
	}

	return matches;
}

// Mostly for testing; doesn't return IDs.
std::vector<size_t> rabin_karp::find_matching_indices(
	const std::vector<char> & haystack) const {

	std::vector<search_result> results = find_matches(haystack);
	std::vector<size_t> match_indices;

	for (const search_result & result: results) {
		match_indices.push_back(result.idx);
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

	rabin_karp rk(str_to_vec(needle), 0);
	if (rk.find_matches(str_to_vec(haystack)).size() > 0) {
		throw std::logic_error("Rabin-Karp test: Found a match "
			"where nowhere was expected");
	}

	// Test matching something.
	needle = "hello";
	rk = rabin_karp(str_to_vec(needle), 0);

	if (!vec_eq(rk.find_matching_indices(str_to_vec(haystack)), {0, 24, 45})) {
		throw std::logic_error("Rabin-Karp test: unexpected output"
			" from three match haystack.");
	}

	// Test needles that are too large to fit in a hash_int.
	needle = "supercalifragilisticexpialidocious";
	haystack = "supercalifragilistic but not "
		"supercalifragilisticexpialidociousexpialidocious";

	rk = rabin_karp(str_to_vec(needle), 0);
	if (!vec_eq(rk.find_matching_indices(str_to_vec(haystack)), {29})) {
		throw std::logic_error("Rabin-Karp test: unexpected output"
			" with very large needle.");
	}

	// Test two needles and IDs.
	std::string needletwo = "listic but not";
	rk = rabin_karp(str_to_vec(needletwo), 1);
	rk.add(str_to_vec(needle), 0);

	if (!vec_eq(rk.find_matches(str_to_vec(haystack)),
		{search_result(14, 1), search_result(29, 0)})) {
		throw std::logic_error("Rabin-Karp test: unexpected output"
			" with two needles.");
	}

	// Test negative values.
	std::vector<char> haystack_nv = {-1, -1, 0, 1, 1, -1, 0, -1},
		needle_nv = {1, 1, -1, 0, -1};
	rk = rabin_karp(needle_nv, 0);

	if (!vec_eq(rk.find_matching_indices(haystack_nv), {3})) {
		throw std::logic_error("Rabin-Karp test: could not find needle"
			" with negative values.");
	}

}