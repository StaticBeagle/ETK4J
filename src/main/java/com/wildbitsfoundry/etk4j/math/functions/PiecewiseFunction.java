package com.wildbitsfoundry.etk4j.math.functions;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.extrapolation.Extrapolator;

public abstract class PiecewiseFunction implements UnivariateFunction {
	private int _numberOfSegments;
	protected double[] _x = null;
	protected double[] coefficients = null;
	// TODO _coefs name doesn't make sense
	Extrapolator _extrapolator;

	private double _x0, _xn;

	protected PiecewiseFunction(double[] x) {
		_x = x;
		_x0 = x[0];
		_numberOfSegments = _x.length - 1;
		_xn = x[_numberOfSegments];
		// TODO check Default extrapolator is throw
		_extrapolator = xi -> {
			throw new IndexOutOfBoundsException(String.format("x is outside of [%.4f, %.4f]", _x0, _xn));
		};
	}

	protected void setExtrapolator(Extrapolator extrapolator) {
		_extrapolator = extrapolator;
	}

	public int findSegmentIndex(double x) {
		int index = Arrays.binarySearch(_x, x);
		return index < 0.0 ? -(index + 2) : Math.min(index, _x.length - 2);
	}

	public int getNumberOfSegments() {
		return _numberOfSegments;
	}

	@Override
	public final double evaluateAt(double x) {
		x += 0.0;	// convert -0.0 to 0.0
		if (x >= _x0 && x <= _xn) {
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
		return this.getSegment(_x.length - 2);
	}

	protected double extrapolate(double x) {
		return _extrapolator.extrapolate(x);
	}

	public double[] getBreaks() {
		return Arrays.copyOf(_x, _x.length);
	}
}
