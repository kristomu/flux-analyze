#pragma once
#include <vector>
#include <cstdio>

// This is used by ordinal search, but also by preamble.h.
// Since ordinal_search refers to preamble, this definition has to be
// put in another file. Kinda hacky, TODO cleanup include spaghetti?

// Ordinal search sequences may actually match one or two
// fluxes into the true search pattern, depending on whether
// the pattern starts with zeroes.

struct ordinal_pattern {
	std::vector<char> needle;
	size_t offset;
};