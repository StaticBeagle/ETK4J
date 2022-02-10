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
	 * Creates a new Nearest Neighbor piecewise interpolant.
	 * @param x The array of abscissas. A copy of this array is made internally.
	 * @param y THe array or ordinates.
	 * @return A newly created Nearest Neighbor piecewise interpolant.
	 */
	public static NearestNeighbor newNearestNeighbor(double[] x, double[] y) {
		return newNearestNeighborInPlace(Arrays.copyOf(x, x.length), y);
	}

	/**
	 * Creates a new Nearest Neighbor piecewise interpolant.
	 * @param x The array of abscissas.
	 * @param y THe array or ordinates.
	 * @return A newly created Nearest Neighbor piecewise interpolant.
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
	 * @return The result of evaluating the interpolant piecewise function.
	 */
	@Override
	public double evaluateAt(int index, double x) {
		double t = (x - this.x[index]) / (this.x[index + 1] - this.x[index]);
		index += t >= 0.5 ? 1 : 0;
		return y[index];
	}

	/**
	 * Get a given segment of the piecewise function.
	 * @param index The index of the segment to fetch.
	 * @return A {@link UnivariateFunction} which represents the underlying segment of the piecewise function.
	 */
	@Override
	public UnivariateFunction getSegment(int index) {
		final double yi = this.evaluateAt(index, x[index]);
		UnivariateFunction fn = x -> yi;
		return fn;
	}

	public static void main(String[] args) {
		double[] x = NumArrays.linSteps(0, 4, 0.5);
		double[] y = Arrays.stream(x).map(v -> v * v).toArray();

		NearestNeighbor nh = NearestNeighbor.newNearestNeighbor(x, y);
		double[] xi = NumArrays.linSteps(0, 4, 0.05);
//		System.out.println();
//		for(double d : xi) {
//			System.out.printf("x: %f, y: %f%n", d, nh.evaluateAt(d));
//		}
		System.out.println(Arrays.toString(nh.evaluateAt(xi)));
		nh.evaluateAt(xi);
		xi = NumArrays.linSteps(0, 10, 0.05);
		x = NumArrays.concatenateAll(NumArrays.linSteps(0, 4.5, 0.5), new double[] {4.99}, NumArrays.linSteps(5, 10, 0.5));// [0:.5:4.5,4.99,5:.5:10];
		y = Arrays.stream(x).map(v -> Math.sin(2.0 * Math.PI * v / 5.0) - (v >= 5 ? 1 : 0)).toArray();
		nh = NearestNeighbor.newNearestNeighbor(x, y);
		//xf = 0:0.05:10;                yf = sin (2*pi*xf/5) - (xf >= 5);
		System.out.println(Arrays.toString(nh.evaluateAt(xi)));

		LinearSpline ls = LinearSpline.newLinearSpline(x, y);
		System.out.println(Arrays.toString(ls.evaluateAt(xi)));
	}
}
