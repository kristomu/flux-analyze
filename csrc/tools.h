#pragma once
#include <vector>
#include <algorithm>

// https://stackoverflow.com/a/14612943
int sign(int x);

template<typename T> double median(std::vector<T> vec) {
	// This will only be called for small vectors, so just sort.
	std::sort(vec.begin(), vec.end());

	if (vec.size() % 2 == 0) {
		return 0.5 * (vec[vec.size()/2-1] + vec[vec.size()/2]);
	} else {
		return vec[vec.size()/2];
	}
}
