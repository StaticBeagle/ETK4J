package com.wildbitsfoundry.etk4j.math.calculus;

import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;

public final class Integrals {
	private Integrals() {
	}

	public static double[] cummulativeTrapz(double... a) {
		final int length = a.length;
		double[] result = new double[length];
		for (int i = 1; i < length; ++i) {
			result[i] = result[i - 1] + (a[i] + a[i - 1]) * 0.5;
		}
		return result;
	}

	public static double trapz(double... a) {
		final int length = a.length;
		double result = 0.0;
		for (int i = 1; i < length; ++i) {
			result += (a[i] + a[i - 1]) * 0.5;
		}
		return result;
	}

	public static double simpson(double... a) {
		final int length = a.length;
		if (length % 2 != 0) {
			throw new IllegalArgumentException("The number of elements must be even");
		}

		double even = 0, odd = 0;
		for (int i = 0; i < length;) {
			odd += a[i++];
			even += a[i++];
		}
		return 1.0 / 3.0 * (a[0] + 2.0 * even + 4.0 * odd + a[length - 1]);
	}

	public static double simpson(UnivariateFunction func, double a, double b, double n) {
		if (n % 2 != 0) {
			throw new IllegalArgumentException("The number of elements must be even");
		}

		double even = 0, odd = 0;
		double h = (b - a) / n;
		for (int i = 1; i < n; i = i + 2) {
			odd += func.evaluateAt(a + i * h);
		}
		for (int i = 2; i < n - 1; i = i + 2) {
			even += func.evaluateAt(a + i * h);
		}
		return h / 3.0 * (func.evaluateAt(a) + 2.0 * even + 4.0 * odd + func.evaluateAt(b));
	}
}
