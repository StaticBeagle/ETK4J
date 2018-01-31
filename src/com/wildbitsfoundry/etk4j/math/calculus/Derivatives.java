package com.wildbitsfoundry.etk4j.math.calculus;

import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;

public final class Derivatives {
	private Derivatives() {
	}

	public static double forwardDifference(UnivariateFunction fn, double x, double h) {
		return (fn.evaluateAt(x + h) - fn.evaluateAt(x)) / h;
	}

	public static double forwardDifference3Points(UnivariateFunction fn, double x, double h) {
		double num = -fn.evaluateAt(x + 2 * h) + 4 * fn.evaluateAt(x + h) - 3 * fn.evaluateAt(x);
		return num / (2 * h);
	}

	public static double forwardDifference5Points(UnivariateFunction fn, double x, double h) {
		double num = -25 * fn.evaluateAt(x) + 48 * fn.evaluateAt(x + h) - 36 * fn.evaluateAt(x + 2 * h)
				+ 16 * fn.evaluateAt(x + 3 * h) - 3 * fn.evaluateAt(x + 4 * h);
		return num / (12 * h);
	}

	public static double backwardDifference(UnivariateFunction fn, double x, double h) {
		return (fn.evaluateAt(x) - fn.evaluateAt(x - h)) / h;
	}

	public static double backwardDifference3Points(UnivariateFunction fn, double x, double h) {
		double num = fn.evaluateAt(x - 2 * h) - 4 * fn.evaluateAt(x - h) + 3 * fn.evaluateAt(x);
		return num / (2 * h);
	}

	public static double backwardDifference5Points(UnivariateFunction fn, double x, double h) {
		double num = 25 * fn.evaluateAt(x) - 48 * fn.evaluateAt(x - h) + 36 * fn.evaluateAt(x - 2 * h)
				- 16 * fn.evaluateAt(x - 3 * h) + 3 * fn.evaluateAt(x - 4 * h);
		return num / (12 * h);
	}

	public static double centeredDifference(UnivariateFunction fn, double x, double h) {
		return (fn.evaluateAt(x + h) - fn.evaluateAt(x - h)) / (2 * h);
	}

	public static double centeredDifference5Points(UnivariateFunction fn, double x, double h) {
		double num = 8 * (fn.evaluateAt(x + h) - fn.evaluateAt(x - h)) - fn.evaluateAt(x + 2 * h)
				+ fn.evaluateAt(x - 2 * h);
		return num / (12 * h);
	}

	public static double centeredDifference7Points(UnivariateFunction fn, double x, double h) {
		double num = 45 * (fn.evaluateAt(x + h) - fn.evaluateAt(x - h))
				+ 9 * (fn.evaluateAt(x - 2 * h) - fn.evaluateAt(x + 2 * h)) + fn.evaluateAt(x + 3 * h)
				- fn.evaluateAt(x - 3 * h);
		return num / (60 * h);
	}
}
