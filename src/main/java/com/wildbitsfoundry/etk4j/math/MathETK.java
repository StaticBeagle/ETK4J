package com.wildbitsfoundry.etk4j.math;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

public final class MathETK {
	private MathETK() {
	}

	/***
	 * Hypotenuse without under/overflow
	 * @param a
	 * @param b
	 * @return sqrt(a^2 + b^2)
	 */
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

	/***
	 * Inverse hyperbolic sine
	 * @param x
	 * @return
	 */
	public static double asinh(double x) {
		return x >= 0 ? Math.log(x + Math.sqrt(x * x + 1.0)) : -Math.log(-x + Math.sqrt(x * x + 1.0));
	}

	public static double acosh(double x) {
		return Math.log(x + Math.sqrt(x * x - 1.0));
	}

	static Complex atanh(double x) {
		double z = Math.abs(x);
		if(z < 1) {	// -1 < x < 1
			double at = 0.5 * Math.log((1 + x) / (1 - x));
			Math.copySign(at, x);
			return Complex.fromReal(at);
		}
		Complex at = Complex.fromImaginary(x).atan();
		at.divideEquals(Complex.fromImaginary(1.0));
		return at;
	}
	
	/***
	 * Round towards zero
	 * @param x
	 * @return the value of x rounded to the nearest integer toward zero
	 */
	public static double fix(double x) {
		return x <= 0 ? Math.ceil(x) : Math.floor(x);
	}
	
	/***
	 * Remainder after division
	 * @param a
	 * @param b
	 * @return
	 */
	public static double rem(double a, double b) {
		if(b == 0) {
			return Double.NaN;
		}
		return a - b * fix(a / b);
	}

	/***
	 * Find the number rounded down to a multiple of the threshold.
	 * @param x
	 * @param threshold
	 * @return
	 */
	public static double floor(double x, double threshold) {
		double result = 0.0d;
		result = x / threshold;
		return Math.floor(result) * threshold;
	}

	/***
	 * Find the number rounded up to a multiple of the threshold.
	 * @param x
	 * @param threshold
	 * @return
	 */
	public static double ceil(double x, double threshold) {
		double result = 0.0d;
		result = x / threshold;
		return Math.ceil(result) * threshold;
	}

	/***
	 *
	 * @param a
	 * @param b
	 * @param absTol
	 * @param relTol
	 * @return {@code Math.abs(a - b) <= absTol + relTol * Math.abs(b)}
	 */
	public static boolean isClose(double a, double b, double absTol, double relTol) {
		return Math.abs(a - b) <= absTol + relTol * Math.abs(b);
	}

	public static boolean isClose(double a, double b, double absTol) {
		return isClose(a, b, absTol, 1e-5);
	}

	public static boolean isClose(double a, double b) {
		return isClose(a, b, 1e-8, 1e-5);
	}

	public static double logb(double x, int base) {
		return Math.log(x) / Math.log(base);
	}

	/***
	 * Round to the next even number.
	 * @param d number to be rounded.
	 * @return the input number rounded to the next even number.
	 */
	public static long roundEven(double d) {
		return Math.round(d / 2) * 2;
	}
}
