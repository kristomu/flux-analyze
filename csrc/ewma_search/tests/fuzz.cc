#include <vector>
#include <stdexcept>

#include <math.h>
#include <stdlib.h>

#include <gtest/gtest.h>

#include "ewma_search/interval.h"
#include "ewma_search/ewma.h"
#include "ewma_search/ewma_search.h"

#include "test.h"

extern "C" int LLVMFuzzerTestOneInput(const char *Data, long long Size) {
	return test_fuzz_input(Data, Size);
}

int main() {
	// Uncomment and add stuff from fuzz tester here when diagnosing
	// such crashes. Remember to uncomment when compiling against libFuzzer.

	test_fuzz_input("\377\377\377\377\377\377\265\275\031\031\030\031\025\000\000\000\025\000\000\000\016\"#\000\004\003\377", 27);
	test_fuzz_input("\377`\"\031I\377\315\274\031\031\030\031\025\000\000\000\025\000\000\000\016\"#\000\004\003\377", 28);
	test_fuzz_input("\377`\"3\000\346\261\274\031\031\030\031\025\000\000\000\025\000\000\000\016\"#\000\004\003\377", 28);
	test_fuzz_input("\377`\"\031I\377\315\274\031\031\030\031\025\000\000\000\025\000\000\000=\"#\000\004\003\377", 27);
	test_fuzz_input("\377`\"\315\377\031\246\274\031\031\030\031\025\000\000\000\025\000\000\000\016\"#\000\004\003\377", 27);
	test_fuzz_input("\377`\"\315\377\031\275\274\031\031\030\031\025\000\000\000\025\000\000\000\016\"#\000\004\003\377", 27);
	test_fuzz_input("\377@I\"\377\031\315\274\031\031\030\031\025\000\000\000\025\000\000\000\016\"#\000\004\003\377", 27);

	std::cout << "OK\n";
}