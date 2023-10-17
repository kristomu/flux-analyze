#pragma once

#include <vector>
#include <cstdint>

// Implement the FNV-1a hash. It's used as a non-cryptographic hash
// to show an alternate hash of data address marks, which should help
// distinguish DAMs with different decodes that happen to have the same
// CRC.

// This code is copied from
// https://gist.github.com/sgsfak/9ba382a0049f6ee885f68621ae86079b

/*
 * The FNV Hash, or more precisely the "FNV-1a alternate algorithm"
 * See: http://www.isthe.com/chongo/tech/comp/fnv/
 *      https://en.wikipedia.org/wiki/Fowler–Noll–Vo_hash_function
 */
uint32_t fnv32_hash(const std::vector<char> & str) {

	/* See the FNV parameters at www.isthe.com/chongo/tech/comp/fnv/#FNV-param */
	const uint32_t FNV_32_PRIME = 0x01000193; /* 16777619 */

	uint32_t h = 0x811c9dc5; /* 2166136261 */

	for (size_t i = 0; i < str.size(); ++i) {
		/* xor the bottom with the current octet */
		h ^= (unsigned char)str[i];
		/* multiply by the 32 bit FNV magic prime mod 2^32 */
		h *= FNV_32_PRIME;
	}

	return h;
}