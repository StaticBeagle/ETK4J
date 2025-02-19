package com.wildbitsfoundry.etk4j.math;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

public final class MathETK {
	private MathETK() {
	}

	/***
	 * Hypotenuse without under/overflow.
	 * @param a The a value
	 * @param b The b value
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
	 * Inverse hyperbolic sine.
	 * @param x The value at where to calculate the function
	 * @return the inverse hyperbolic sine
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
	 * Round towards zero.
	 * @param x The value at where to calculate the function
	 * @return the value of x rounded to the nearest integer toward zero
	 */
	public static double fix(double x) {
		return x <= 0 ? Math.ceil(x) : Math.floor(x);
	}
	
	/***
	 * Remainder after division.
	 * @param a The dividend
	 * @param b The divisor
	 * @return The remainder of a ÷ b
	 */
	public static double rem(double a, double b) {
		if(b == 0) {
			return Double.NaN;
		}
		return a - b * fix(a / b);
	}

	/***
	 * Find the number rounded down to a multiple of the threshold.
	 * @param x The input value
	 * @param threshold The threshold value
	 * @return The value of x rounded down to the nearest multiple of the threshold
 	 */
	public static double floor(double x, double threshold) {
		double result = 0.0d;
		result = x / threshold;
		return Math.floor(result) * threshold;
	}

	/***
	 * Find the number rounded up to a multiple of the threshold.
	 * @param x The input value
	 * @param threshold The threshold value
	 * @return The value of x rounded up to the nearest multiple of the threshold
	 */
	public static double ceil(double x, double threshold) {
		double result = 0.0d;
		result = x / threshold;
		return Math.ceil(result) * threshold;
	}

	/***
	 * Determines whether a number is close to another number within a certain tolerance.
	 * @param a Argument in which to evaluate the function at.
	 * @param b Argument in which to evaluate the function at.
	 * @param absTol The absolute tolerance.
	 * @param relTol The relative tolerance.
	 * @return {@code Math.abs(a - b) <= absTol + relTol * Math.abs(b)}
	 * @see <a href="https://numpy.org/doc/stable/reference/generated/numpy.isclose.html">isClose</a>
	 */
	public static boolean isClose(double a, double b, double absTol, double relTol) {
		return Math.abs(a - b) <= absTol + relTol * Math.abs(b);
	}

	/***
	 * Determines whether a number is close to another number within a certain tolerance.
	 * @param a Argument in which to evaluate the function at.
	 * @param b Argument in which to evaluate the function at.
	 * @param absTol The absolute tolerance.
	 * @param relTol The relative tolerance.
	 * @return {@code Math.abs(a - b) <= absTol + relTol * Math.abs(b)}
	 * @see <a href="https://numpy.org/doc/stable/reference/generated/numpy.isclose.html">isClose</a>
	 */
	public static boolean isClose(Complex a, Complex b, double absTol, double relTol) {
		return a.subtract(b).abs() <= absTol + relTol * b.abs();
	}

	/***
	 * Determines whether a number is close to another number within a certain tolerance. <br>
	 * This calls {@link #isClose(double, double, double, double) with relative tolerance of 1e-5}.
	 * @param a Argument in which to evaluate the function at.
	 * @param b Argument in which to evaluate the function at.
	 * @param absTol The absolute tolerance.
	 * @return {@code Math.abs(a - b) <= absTol + 1e-5 * Math.abs(b)}
	 * @see <a href="https://numpy.org/doc/stable/reference/generated/numpy.isclose.html">isClose</a>
	 */
	public static boolean isClose(double a, double b, double absTol) {
		return isClose(a, b, absTol, 1e-5);
	}

	/***
	 * Determines whether a number is close to another number within a certain tolerance. <br>
	 * This calls {@link #isClose(double, double, double, double) with absolute tolerance of 1e-8 and <br>
	 * relative tolerance of 1e-5}.
	 * @param a Argument in which to evaluate the function at.
	 * @param b Argument in which to evaluate the function at.
	 * @return {@code Math.abs(a - b) <= 1e-8 + 1e-5 * Math.abs(b)}
	 * @see <a href="https://numpy.org/doc/stable/reference/generated/numpy.isclose.html">isClose</a>
	 */
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

	public static class FRexpResult {
		public int exponent = 0;
		public double mantissa = 0.0;
	}

	/**
	 * Find the mantissa and exponent of a number.
	 *
	 * @param value Breaks the floating point number value into its binary significand and an integral exponent for 2.
	 * @return The mantissa and exponent of a number such that number = m * 2^e.
	 * @see <a href="https://stackoverflow.com/a/3946294/6383857">https://stackoverflow.com/a/3946294/6383857</a>
	 */
	public static FRexpResult frExp(double value) {
		FRexpResult result = new FRexpResult();

		result.exponent = 0;
		result.mantissa = 0;

		if (value == 0.0) {
			return result;
		}
		if (Double.isNaN(value)) {
			result.mantissa = Double.NaN;
			result.exponent = -1;
			return result;
		}
		if (Double.isInfinite(value)) {
			result.mantissa = value;
			result.exponent = -1;
			return result;
		}

		long bits = Double.doubleToLongBits(value);
		double realMantissa = 1.0;

		boolean neg = (bits < 0);
		int exponent = (int) ((bits >> 52) & 0x7ffL);
		long mantissa = bits & 0xfffffffffffffL;

		if (exponent == 0) {
			exponent++;
		} else {
			mantissa = mantissa | (1L << 52);
		}

		// bias the exponent - actually biased by 1023.
		// we are treating the mantissa as m.0 instead of 0.m
		//  so subtract another 52.
		exponent -= 1075;
		realMantissa = mantissa;

		// normalize
		while (realMantissa > 1.0) {
			mantissa >>= 1;
			realMantissa /= 2.0;
			exponent++;
		}

		if (neg) {
			realMantissa = realMantissa * -1;
		}

		result.exponent = exponent;
		result.mantissa = realMantissa;

		return result;
	}

	/**
	 * Combinations of picking k unordered outcomes from n possibilities
	 * @param n the number of things
	 * @param k the number of elements taken
	 * @return The total number of combinations
	 */
	public static long combinations(long n, long k) {
		return factorial(n) / (factorial(k) * factorial(n - k));
	}

	/**
	 * Factorial of a number
	 * @param n the number
	 * @return {@code n!}
	 */
	public static long factorial(long n) {
		if(n == 0) {
			return 1;
		}
		return n * factorial(n - 1);
	}
}
