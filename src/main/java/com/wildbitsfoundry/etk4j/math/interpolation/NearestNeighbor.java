package com.wildbitsfoundry.etk4j.math.interpolation;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.extrapolation.Extrapolators;
import com.wildbitsfoundry.etk4j.math.functions.PiecewiseFunction;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import com.wildbitsfoundry.etk4j.util.NumArrays;

import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkMinXLength;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkXYDimensions;

/**
 * The {@code NearestNeighbor} class implements 1d interpolation using the Nearest Neighbor method.
 * @see <a href="https://en.wikipedia.org/wiki/Nearest-neighbor_interpolation">Nearest Neighbor interpolation</a>
 */
public class NearestNeighbor extends PiecewiseFunction {
	
	private double[] y;

	protected NearestNeighbor(double[] x, double[] y) {
		super(x);
		this.y = Arrays.copyOf(y, y.length);

		double x0 = x[0];
		double xn = x[x.length - 1];
		double y0 = this.evaluateAt(x0);
		double yn = this.evaluateAt(xn);
		this.setExtrapolator(new Extrapolators.ClampToEndPointExtrapolator(x0, xn, y0, yn));
	}

	/**
	 * Creates a new Nearest Neighbor instance.
	 * @param x The array of x coordinates. The values in this array must be unique and strictly increasing.
	 *          A copy of this array is made internally.
	 * @param y THe array or y coordinates.
	 * @return A newly created Nearest Neighbor interpolant.
	 */
	public static NearestNeighbor newNearestNeighbor(double[] x, double[] y) {
		return newNearestNeighborInPlace(Arrays.copyOf(x, x.length), y);
	}

	/**
	 * Creates a new Nearest Neighbor instance.
	 * @param x The array of abscissas.
	 * @param y THe array or ordinates.
	 * @return A newly created Nearest Neighbor interpolant.
	 */
	public static NearestNeighbor newNearestNeighborInPlace(double[] x, double[] y) {
		checkXYDimensions(x, y);
		checkMinXLength(x, 2);
		NearestNeighbor nh = new NearestNeighbor(x, y);
		return nh;
	}

	/**
	 * Evaluate the piecewise function.
	 * @param index The index of the segment to evaluate the function at.
	 * @param x The value at which to evaluate the function.
	 * @return The result of evaluating the piecewise function.
	 */
	@Override
	public double evaluateAt(int index, double x) {
		double t = (x - this.x[index]) / (this.x[index + 1] - this.x[index]);
		index += t >= 0.5 ? 1 : 0;
		return y[index];
	}

	/**
	 * Get a given segment of the piecewise function.
	 * @param index The index of the segment to get.
	 * @return A {@link UnivariateFunction} which represents the underlying segment of the piecewise function.
	 */
	@Override
	public UnivariateFunction getSegment(int index) {
		final double yi = this.evaluateAt(index, x[index]);
		UnivariateFunction fn = x -> yi;
		return fn;
	}
}
