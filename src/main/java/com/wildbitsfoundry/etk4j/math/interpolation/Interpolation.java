package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.math.curvefitting.CurveFitting;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkXYDimensions;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkMinXLength;

/**
 * The {@code Interpolation} class contains methods to perform interpolation at a given argument.
 */
public final class Interpolation {

	private Interpolation() {
	}

	/**
	 * Linear interpolation. Creates a line and evaluates that line at the value specified by {@code xi}.
	 * @param x0 Lower x coordinate.
	 * @param x1 Upper x coordinate.
	 * @param y0 Lower ordinate.
	 * @param y1 Upper ordinate.
	 * @param xi Argument at which to evaluate the function.
	 * @return {@code y(xi) = m * xi + b}.
	 */
	public static double linear(double x0, double x1, double y0, double y1, double xi) {
		if(xi < x0 || xi > x1) {
			throw new IllegalArgumentException("xi must be within [x0, x1].");
		}
		double hx = x1 - x0;
		double t = (xi - x0) / hx;
		return (y1 - y0) * t + y0;
	}
	// TODO check all the See in the project and change to @see if needed
	/**
	 * Neville's interpolation algorithm. This array must be monotonically increasing.
	 * @param x Array of x coordinates.
	 * @param y Array of ordinates.
	 * @param xi Argument at which to evaluate the function.
	 * @return The result of Neville's algorithm evaluated at {@code xi}.
	 * @see <a href="https://mathworld.wolfram.com/NevillesAlgorithm.html">Neville's algorithm</a>
	 */
	public static double neville(double[] x, double[] y, double xi) {
		if(xi < x[0] || xi > x[x.length - 1]) {
			throw new IllegalArgumentException(String.format("xi must lie within the bounds of x [%.4f, %.4f].",
					x[0], x[x.length - 1]));
		}
		checkXYDimensions(x, y);
		int length = x.length;
		double[][] N = new double[length][length];

		// Initializing first column
		for (int i = 0; i < length; ++i) {
			N[i][0] = y[i];
		}
		// Neville's method.
		for (int i = 1; i < length; ++i) {
			for (int j = 1; j <= i; ++j) {
				N[i][j] = ((xi - x[i - j]) * (N[i][j - 1]) - (xi - x[i]) * (N[i - 1][j - 1])) / (x[i] - x[i - j]);
			}
		}
		return N[length - 1][length - 1];
	}

	/**
	 * Quadratic (Parabola) interpolation.
	 * @param x0 Left x coordinate.
	 * @param x1 Middle x coordinate.
	 * @param x2 Right x coordinate.
	 * @param y0 Left y coordinate.
	 * @param y1 Middle y coordinate.
	 * @param y2 Right y coordinate.
	 * @param xi Argument at which to evaluate the function.
	 * @return Performs a parabolic fit and evaluates the function at the value specified by {@code xi}.
	 */
	public static double quadratic(double x0, double x1, double x2, double y0, double y1, double y2, double xi) {
		if(xi < x0 || xi > x2) {
			throw new IllegalArgumentException("xi must be within [x0, x2].");
		}
		double[] parabola = CurveFitting.parabola(x0, x1, x2, y0, y1, y2);
		return parabola[2] + xi * (parabola[1] + xi * parabola[0]);
	}

	/**
	 * Spline interpolation.
	 * @param x Array of x coordinates. This array must be monotonically increasing.The values in this array must be
	 *          unique and strictly increasing.
	 * @param y Array of y coordinates.
	 * @param xi Argument at which to evaluate the spline.
	 * @return If the number of points in x ie equal to 2 it performs linear interpolation. If the number of points
	 * in x is equal to 3 it performs quadratic interpolation. If the number of points is greater than 3 it creates
	 * a cubic spline with {@link CubicSpline#newNotAKnotSpline(double[], double[])} with not-a-knot conditions and
	 * evaluates given spline at {@code xi}.
	 */
	public static double spline(double[] x, double[] y, double xi) {
		checkXYDimensions(x, y);
		checkMinXLength(x, 2);
		
		final int length = x.length;
		if (length == 2) {
			return linear(x[0], x[1], y[0], y[1], xi);
		}
		if(length == 3) {
			return quadratic(x[0], x[1], x[2], y[0], y[1], y[2], xi);
		}
		return CubicSpline.newCubicSpline(x, y).evaluateAt(xi);
	}
}
