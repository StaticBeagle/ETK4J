package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.math.functions.DifferentiableFunction;
import com.wildbitsfoundry.etk4j.math.functions.IntegrableFunction;
import com.wildbitsfoundry.etk4j.math.functions.PiecewiseFunction;

public abstract class Spline extends PiecewiseFunction implements DifferentiableFunction, IntegrableFunction {

	protected Spline(double[] x, double y0, double yn) {
		super(x, y0, yn);
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
