// Implementation of the Rabin-Karp search algorithm. I might
// add an error-tolerant version later (the principle is easy;
// just generate hashes for each version with errors).

// Worst case complexity is O(nm). Aho-Corasick would be better
// and not restricted to a given string length, but also a much
// greater pain to write.

#include <vector>
#include <unordered_map>

#include <stdlib.h>

typedef int64_t hash_int;

class rabin_karp {
	private:
		std::unordered_map<hash_int,
			std::vector<char> > needles_by_hash;

		// These are magic constants; they should be
		// reasonably good. The requirements are that
		// both the factor and modulus should be prime,
		// that the modulus is large, and that factor *
		// modulus * max haystack value (here 255) is
		// less than the maximum of the type.
		static const hash_int factor = 3;
		static const hash_int modulus = 9007199254740997;

		size_t min_needle_length;
		hash_int leading_char_eliminator;

		bool brute_force_search(std::vector<char>::const_iterator pos,
			const std::vector<char> & needle) const;

	public:
		rabin_karp(size_t min_needle_length_in);
		rabin_karp(const std::vector<char> & needle);
		void add(const std::vector<char> & needle);
		std::vector<size_t> find_matches(
			const std::vector<char> & haystack) const;
};

// Unit test. I should probably get a proper test framework for this...
void test_rabin_karp();