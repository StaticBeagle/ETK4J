package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.math.functions.DifferentiableFunction;
import com.wildbitsfoundry.etk4j.math.functions.IntegrableFunction;
import com.wildbitsfoundry.etk4j.math.functions.PiecewiseFunction;

public abstract class Spline extends PiecewiseFunction implements DifferentiableFunction, IntegrableFunction {

	public enum SplineType {
		CUBIC,
		LINEAR,
		MONOTONIC
	}
	
	private double[] _indefiniteIntegral = null;
	protected double[] _coefs = null;
	private final int _order;

	protected Spline(double[] x, int order, double y0, double yn) {
		super(x, y0, yn);
		_order = order;
	}
	
	protected abstract double evaluateDerivativeAt(int i, double t);

	protected abstract double evaluateAntiDerivativeAt(int i, double t);

	public int getOrder() {
		return _order;
	}
	
	@Override
	protected double getValueAt(int i, double x) {
		
		double t = x - _x[i];
		i *= _order;
		double result = 0;
		for (int j = 0; j < _order; ++j) {
			result = result * t + _coefs[i++];
		}
		return result;
	}
	
	@Override
	public double differentiate(double x) {
		int i = this.findLeftIndex(x);
		double t = x - _x[i];
		return this.evaluateDerivativeAt(i, t);
	}

	@Override
	public double integrate(double a, double b) {
		if (a < _x[0] || b > _x[_x.length - 1]) {
			throw new IllegalArgumentException(
					String.format("The spline is not defined outside of [%.4g, %.4g]", _x[0], _x[_x.length - 1]));
		}
		return integrate(b) - integrate(a);
	}

	private double integrate(double x) {
		// Lazy creating of values
		if (_indefiniteIntegral == null) {
			this.calcuateIntegral();
		}

		int i = this.findLeftIndex(x);
		double t = x - _x[i];
		return _indefiniteIntegral[i] + this.evaluateAntiDerivativeAt(i, t);
	}

	private void calcuateIntegral() {
		// Lazy creating the values
		if (_indefiniteIntegral != null) {
			return;
		}

		final int size = _coefs.length / _order;
		_indefiniteIntegral = new double[size];
		for (int i = 0; i < size - 1; ++i) {
			double t = _x[i + 1] - _x[i];
			_indefiniteIntegral[i + 1] = _indefiniteIntegral[i] + this.evaluateAntiDerivativeAt(i, t);
		}
	}
	
	public Spline newSpline(double[] x, double[] y, SplineType type) {
		switch (type) {
		case CUBIC:
			return CubicSpline.newCubicSpline(x, y);
		case LINEAR:
			return LinearSpline.newLinearSpline(x, y);
		case MONOTONIC:
			return CubicSpline.newAkimaSpline(x, y);
		default:
			throw new IllegalStateException();
		}
	}
	
	protected static void checkXYDimensions(double[] x, double[] y) {
		if (x.length != y.length) {
			throw new IllegalArgumentException("x and y dimensions must be the same");
		}
	}

	protected static void checkMinkXLength(double[] x, int size) {
		if (x.length < size) {
			throw new IllegalArgumentException(String.format("x length must be >= %d", size));
		}
	}
}
