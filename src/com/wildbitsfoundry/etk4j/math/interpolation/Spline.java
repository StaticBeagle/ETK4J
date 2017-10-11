package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.math.functions.DifferentiableFunction;
import com.wildbitsfoundry.etk4j.math.functions.IntegrableFunction;
import com.wildbitsfoundry.etk4j.math.functions.PiecewiseFunction;

public abstract class Spline extends PiecewiseFunction implements DifferentiableFunction, IntegrableFunction {

	private double[] _indefiniteIntegral = null;
	protected double[] _coefs = null;
	
	protected Spline(double[] x, double y0, double yn) {
		super(x, y0, yn);
	}
	
	protected abstract double evaluateAntiDerivativeAt(int index, double t);
	
	public abstract int getOrder();
	
	@Override
	public double integrate(double a, double b) {
		return integrate(b) - integrate(a);
	}
	
	private double integrate(double x) {
		// Lazy creating of values
		if(_indefiniteIntegral == null) {
			this.calcuateIntegral();
		}
		
		int i = this.findIndex(x);
		double t = x - _x[i];
		return _indefiniteIntegral[i] + this.evaluateAntiDerivativeAt(i, t);
	}
	
	
	private void calcuateIntegral() {
		// Lazy creating the values
		if(_indefiniteIntegral != null) {
			return;
		}
		
		final int size = _coefs.length / (this.getOrder() + 1);
		_indefiniteIntegral = new double[size];
		for(int i = 0; i < size - 1; ++i) {
			double t = _x[i + 1] -_x[i];
			_indefiniteIntegral[i + 1] = _indefiniteIntegral[i] + this.evaluateAntiDerivativeAt(i, t);
		}
	}

	protected static void checkXYDimensions(double[] x, double[] y) {
		if(x.length != y.length) {
			throw new IllegalArgumentException("x and y dimensions must be the same");
		}
	}
	
	protected static void checkMinkXLength(double[] x, int size) {
		if(x.length < size) {
			throw new IllegalArgumentException(String.format("x length must be >= %d", size));
		}
	}
}
