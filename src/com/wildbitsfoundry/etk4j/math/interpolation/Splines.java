package com.wildbitsfoundry.etk4j.math.interpolation;

public final class Splines {
	private Splines() {}
	
	public enum SplineType {
		CUBIC,
		LINEAR,
		MONOTONIC
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
