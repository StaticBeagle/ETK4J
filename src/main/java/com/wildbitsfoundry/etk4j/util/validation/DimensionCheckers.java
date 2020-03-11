package com.wildbitsfoundry.etk4j.util.validation;

public final class DimensionCheckers {
	private DimensionCheckers() {
		
	}
	
	public static void checkXYDimensions(double[] x, double[] y) {
		if (x.length != y.length) {
			throw new IllegalArgumentException("x and y dimensions must be the same");
		}
	}

	public static void checkMinXLength(double[] x, int minLength) {
		if (x.length < minLength) {
			throw new IllegalArgumentException(String.format("x length must be >= %d", minLength));
		}
	}
}
