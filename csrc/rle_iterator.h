#pragma once

#include <vector>
#include <iostream>
#include <algorithm>

// This run-length iterator turns a vector like {3, 4, 2} into
// 0 0 0 1 0 0 0 0 1 0 0 1. This is used for the MFM decoder
// because it's much easier to decode to MFM at the transition/no
// transition level than at the run length level.

// The way this works is that the iterator keeps a digit and an
// underlying position. When advancing the iterator, it either
// decrements the digit or advances to the new underlying position
// and sets the digit to the value at that position. Thus if the
// first position is 10, it'll decrement 10 times before advancing,
// hence producing a run of then zeroes.

// TODO: Deal with boundary conditions, either here or in sector_data.cc
// E.g. if our base vector's starting point says "3", then that probably
// means 1 0 0 0 1 ..., not "just" 0 0 0 1. But making the iterator output
// 1 first would make things *very* messy, because it means that an RLE
// iterator has to start "before the start". It's probably better to make
// Rabin-Karp return a hit just after the first A1/C2, i.e. if the sequence is
//  0 1 0 0 0 1 0 0 1 0 0 0 1 0 0 1 with run length values
// ???? ---3--- -2--- ---3--- --2--
// then it should search for "3 2 3 2", create an iterator at this point,
// and then rewind the iterator two steps.

// This means that I need to make an operator--. Do that later.

class run_length_const_iterator {
	using iterator_category = std::input_iterator_tag;
	using value_type        = unsigned int;
	using pointer           = unsigned int*;
	using reference         = unsigned int&;

	private:
		const std::vector<unsigned int>::const_iterator begin, end;
		std::vector<unsigned int>::const_iterator underlying_pos;
		unsigned int digit;

	public:
		run_length_const_iterator(const std::vector<unsigned int> & base_vector)
			: begin(base_vector.begin()), end(base_vector.end()),
			  underlying_pos(begin), digit(base_vector[0]) {}

		run_length_const_iterator(
				const std::vector<unsigned int>::const_iterator & begin_in,
				const std::vector<unsigned int>::const_iterator & end_in,
				const std::vector<unsigned int>::const_iterator pos)
			: begin(begin_in), end(end_in),
			  underlying_pos(pos), digit(0) {

			// Special-case things: if the user is requesting we
			// start at the end of the vector, then the digit must
			// be zero for comparisons to work properly. Otherwise,
			// the digit must be set to whatever's at the position.
			if (pos != end) {
				digit = *pos;
			}
		}

		run_length_const_iterator(
				const std::vector<unsigned int> & base_vector,
				const std::vector<unsigned int>::const_iterator pos)
			: run_length_const_iterator(base_vector.begin(),
				base_vector.end(), pos) {}

		int operator*() const {
			if (underlying_pos == end) {
				throw std::out_of_range("run_length_const_iterator:"
					" tried to access end!");
			}
			if (digit > 0) {
				return 0;
			} else {
				return 1;
			}
		}

		run_length_const_iterator & operator++() {
			// If the digit is nonzero, decrement it by one.
			if (digit != 0) {
				--digit;
				return *this;
			}

			// Otherwise, we need to go to the next element of the
			// base vector and set the digit to its value. Don't do
			// that if we're already at the end, though, because we
			// can go no further.
			if (underlying_pos == end) {
				return *this;
			}

			++underlying_pos;
			digit = *underlying_pos;
			return *this;
		}

		// Postfix increment should return a copy
		run_length_const_iterator operator++(int) {
			run_length_const_iterator tmp = *this;
			++tmp;
			return tmp;
		}

		bool operator!=(const run_length_const_iterator & other) const {
			return other.underlying_pos != underlying_pos ||
				   other.digit != digit;
		}

		bool operator==(const run_length_const_iterator & other) const {
			return !(*this != other);
		}
};