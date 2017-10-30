package com.wildbitsfoundry.etk4j.math;

public final class MathETK {
	private MathETK() {
	}

	public static double hypot(double a, double b) {
		double r;
		if (Math.abs(a) > Math.abs(b)) {
			r = b / a;
			r = Math.abs(a) * Math.sqrt(1 + r * r);
		} else if (b != 0) {
			r = a / b;
			r = Math.abs(b) * Math.sqrt(1 + r * r);
		} else {
			r = 0.0;
		}
		return r;
	}

	public static double asinh(double x) {
		return Math.log(x + Math.sqrt(x * x + 1.0));
	}

	public static double acosh(double x) {
		return Math.log(x + Math.sqrt(x * x - 1.0));
	}

	static double atanh(double x) {
		return 0.5 * Math.log((x + 1.0) / (x - 1.0));
	}
}
