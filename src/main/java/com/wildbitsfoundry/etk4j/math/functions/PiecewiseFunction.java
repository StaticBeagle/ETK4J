package com.wildbitsfoundry.etk4j.math.functions;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.extrapolation.Extrapolator;
import com.wildbitsfoundry.etk4j.math.extrapolation.Extrapolators;

/**
 * The {@code PiecewiseFunction} class describes the mathematical equivalent of a mathematical piecewise function. In
 * other words, a function comprised of segments which are functions themselves.
 */
public abstract class PiecewiseFunction implements UnivariateFunction {
	private int numberOfSegments;
	protected double[] x;
	protected double[] coefficients = null;
	Extrapolator extrapolators;

	private double x0, xn;

	protected PiecewiseFunction(double[] x) {
		this.x = x;
		x0 = x[0];
		numberOfSegments = this.x.length - 1;
		xn = x[numberOfSegments];
		// TODO check Default extrapolator is throw
		extrapolators = new Extrapolators.ThrowExtrapolator(x0, xn);
	}

	protected void setExtrapolator(Extrapolator extrapolator) {
		extrapolators = extrapolator;
	}

	public int findSegmentIndex(double x) {
		int index = Arrays.binarySearch(this.x, x);
		return index < 0.0 ? -(index + 2) : Math.min(index, this.x.length - 2);
	}

	public int getNumberOfSegments() {
		return numberOfSegments;
	}

	@Override
	public final double evaluateAt(double x) {
		x += 0.0;	// convert -0.0 to 0.0
		if (x >= x0 && x <= xn) {
			int index = this.findSegmentIndex(x);
			return this.evaluateAt(index, x);
		}
		return this.extrapolate(x);
	}

	public final double[] evaluateAt(double[] x) {
		final int n = x.length;
		double[] yi = new double[n];
		for(int i = 0; i < n; ++i) {
			yi[i] = this.evaluateAt(x[i]);
		}
		return yi;
	}

	public abstract double evaluateAt(int index, double x);

	public abstract UnivariateFunction getSegment(int index);

	public UnivariateFunction getFirstSegment() {
		return this.getSegment(0);
	}

	public UnivariateFunction getLastSegment() {
		return this.getSegment(x.length - 2);
	}

	protected double extrapolate(double x) {
		return extrapolators.extrapolate(x);
	}

	public double[] getBreaks() {
		return Arrays.copyOf(x, x.length);
	}
}
