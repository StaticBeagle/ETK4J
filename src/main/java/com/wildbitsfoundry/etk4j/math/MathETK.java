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
			return Complex.of(at, 0.0);
		}
		Complex at = Complex.of(0.0, x).atan();
		at.divideEquals(Complex.of(0.0, 1.0));
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
}
