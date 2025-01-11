package com.wildbitsfoundry.etk4j.util.validation;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

public final class DimensionCheckers {
	private DimensionCheckers() {
		
	}
	
	public static void checkXYDimensions(double[] x, double[] y) {
		if (x.length != y.length) {
			throw new IllegalArgumentException("The length of both arrays must match");
		}
	}

	public static void checkXYDimensions(Complex[] x, Complex[] y) {
		if (x.length != y.length) {
			throw new IllegalArgumentException("The length of both arrays must match");
		}
	}

	public static void checkMinXLength(double[] x, int minLength) {
		if (x.length < minLength) {
			throw new IllegalArgumentException(String.format("The length of the array must be >= %d", minLength));
		}
	}
}
