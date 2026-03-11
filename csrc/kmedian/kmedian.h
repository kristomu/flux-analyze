#include <algorithm>
#include <numeric>
#include <math.h>
#include <limits.h>
#include <vector>
#include <iostream>
#include <iterator>
#include <functional>

/* Implement a k * n log n method for optimal 1D K-median.
 * The K-median problem consists of finding k clusters so that the
 * sum of absolute distances from each of the n points to the closest
 * cluster is minimized.
 *
 * GRØNLUND, Allan, et al. Fast exact k-means, k-medians and Bregman
 * divergence clustering in 1d. arXiv preprint arXiv:1701.07204, 2017.
 */

class KCluster
{
public:
	int numPoints;                    /* number of points */
	std::vector<float> points;        /* sorted list of points */
	std::vector<float> cumulativeSum; /* cumulative sum of points */

	KCluster(std::vector<float> points_in);

	/* Returns the median point from the sorted list of points in the supplied
	 * half-open range. */

	float medianPoint(int low, int high) const;

	 /* Grønlund et al.'s CC function for the K-median problem.
	  *
	  * CC(i,j) is the cost of grouping points_i...points_j into one cluster with the
	  * optimal cluster point (the median point). By programming convention, the
	  * interval is half-open and indexed from 0, unlike the paper's convention
	  * of closed intervals indexed from 1.
	  */

	float CC(int i, int j) const;

	/* D_previous is the D vector for (i-1) clusters, or empty if i < 2.
	 * It's possible to do this even faster (and more incomprehensibly).
	 * See Grønlund et al. for that. */

	/* Calculate C_i[m][j] given D_previous = D[i-1]. p. 4 */

	float C_i(int i, const std::vector<float>& D_previous, int m,
			int j) const;

	/* Find the (location of the) minimum value of each row of an nxn matrix in
	 * n log n time, given that the matrix is monotone (by the definition of
	 * Grønlund et al.) and defined by a function M that takes row and column
	 * as parameters. All boundary values are closed, [minRow, maxRow]
	 * p. 197 of
	 * AGGARWAL, Alok, et al. Geometric applications of a matrix-searching
	 * algorithm. Algorithmica, 1987, 2.1-4: 195-208.
	 */

	void monotoneMatrixIndices(std::function<float(int, int)>& M,
			int minRow, int maxRow, int minCol, int maxCol,
			std::vector<int>& Tout,
			std::vector<float>& Dout) const;

	std::vector<float> optimalKMedian(int numClusters) const;
};

/* TEST
 */

bool test();
