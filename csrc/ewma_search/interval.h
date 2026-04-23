#pragma once

// Define a closed interval of values (here, admissible clock
// frequencies).

// TODO? Use infinities?

#include <cstdlib>			// size_t

class interval {
	public:
		double lower_clock = -1, upper_clock = -1;
		bool is_set = false;

		bool valid() {
			return lower_clock <= upper_clock;
		}

		void shrink_bounds(double new_lower_clock,
			double new_upper_clock) {

			if (is_set) {
				lower_clock = std::max(lower_clock,
					new_lower_clock);
				upper_clock = std::min(upper_clock,
					new_upper_clock);
			} else {
				lower_clock = new_lower_clock;
				upper_clock = new_upper_clock;
			}

			is_set = true;
		}

		void grow_bounds(double new_lower_clock,
			double new_upper_clock) {

			if (is_set) {
				lower_clock = std::min(lower_clock,
					new_lower_clock);
				upper_clock = std::max(upper_clock,
					new_upper_clock);
			} else {
				lower_clock = new_lower_clock;
				upper_clock = new_upper_clock;
			}

			is_set = true;
		}

		interval() {}

		interval(double lower_in, double upper_in) {
			shrink_bounds(lower_in, upper_in);
		}

		bool operator==(const interval & other) const {
			return (is_set == other.is_set &&
				lower_clock == other.lower_clock &&
				upper_clock == other.upper_clock);
		}

};
