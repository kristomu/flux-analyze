#include <vector>

// For testing - to determine if the sophisticated search
// algorithm matches iff there is a match. (Up to some margin
// of confidence; test can never decide the matter with certainty,
// of course.)
std::pair<int, interval> brute_force_ewma_search(
	const std::vector<int> & haystack,
	const std::vector<int> & half_clock_needle,
	double alpha);

// A more sophisticated (but slower) version.
std::pair<int, interval> half_brute_force_ewma_search(
	const std::vector<int> & haystack,
	const std::vector<int> & half_clock_needle,
	double alpha);

// Used for fuzz/property testing. If the brute-force
// searches don't match the EWMA search, then at least
// one of them is wrong. This halts (abort()) on mismatch,
// so a fuzz tester can record the error as a crash.
int test_fuzz_input(const char *Data, long long Size);