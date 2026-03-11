#pragma once

#include <functional>
#include <vector>

class KCluster {
	private:
		struct CDF
		{
			unsigned int pointsSoFar = 0;
			unsigned long long sumToHere = 0;
			unsigned int thisValue = 0;
		};

		bool initialized;

		size_t uniquePoints;         /* number of unique points */
		std::vector<struct CDF> cdf; /* cumulative sum of points */

		/* Determines the median point of all points between the ith category
		 * in cdf (inclusive), and the jth (exclusive). It's used to
		 * reconstruct clusters later.
		 */

		double medianPoint(size_t i, size_t j) const;

		/* Find the (location of the) minimum value of each row of an nxn matrix in
		 * n log n time, given that the matrix is monotone (by the definition of
		 * Grønlund et al.) and defined by a function M that takes row and column
		 * as parameters. All boundary values are closed, [minRow, maxRow]
		 * p. 197 of
		 * AGGARWAL, Alok, et al. Geometric applications of a matrix-searching
		 * algorithm. Algorithmica, 1987, 2.1-4: 195-208.
		 */

		void monotoneMatrixIndices(std::function<double(size_t, size_t)>& M,
		        size_t minRow, size_t maxRow, size_t minCol, size_t maxCol,
		        std::vector<size_t>& Tout,
		        std::vector<double>& Dout) const;

		/* If item is the first cdf value with count at
		 * or above i, returns the cumulative sum up to the ith entry of the
		 * underlying sorted points list.
		 */

		unsigned long long cumulativeAt(
			std::vector<struct CDF>::const_iterator item,
			int i) const;

		/* Grønlund et al.'s CC function for the K-median problem.
		 *
		 * CC(i,j) is the cost of grouping points_i...points_j into one cluster
		 * with the optimal cluster point (the median point). By programming
		 * convention, the interval is half-open and indexed from 0, unlike the
		 * paper's convention of closed intervals indexed from 1.
		 *
		 * Note: i and j are indices onto weighted_cdf. So e.g. if i = 0, j = 2 and
		 * weighted cdf is [[2, 0, 0], [4, 2, 1], [5, 2, 2]], then that is the cost
		 * of clustering all points between the one described by the zeroth weighted
		 * cdf entry, and up to (but not including) the last. In other words, it's
		 * CC([0, 0, 1, 1], 0, 5).
		 */

		long double CC(size_t i, size_t j) const;

		/* D_previous is the D vector for (i-1) clusters, or empty if i < 2.
		 * It's possible to do this even faster (and more incomprehensibly).
		 * See Grønlund et al. for that. */

		/* Calculate C_i[m][j] given D_previous = D[i-1]. p. 4 */

		long double C_i(int i, const std::vector<double> & D_previous,
			size_t m, size_t j) const;

		/* Calculates the optimal cluster centres for the points. */

	public:

		void setPoints(const std::vector<int> & inputPoints);

		KCluster(const std::vector<int> & inputPoints) {
			setPoints(inputPoints);
		}

		KCluster() {
			initialized = false;
		}

		std::vector<double> optimalKMedian(int numClusters);

		std::vector<double> optimalKMedian(
			const std::vector<int> & inputPoints,
			int numClusters) {

			setPoints(inputPoints);
			return optimalKMedian(numClusters);
		}

};
