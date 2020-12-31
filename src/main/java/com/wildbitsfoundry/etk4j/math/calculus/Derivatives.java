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

	// http://www.ibiblio.org/technicalc/tiplist/en/files/pdf/tips/tip6_26.pdf
	public static double ridders(UnivariateFunction fn, double x, double h) {
		if(h == 0) {
			throw new IllegalArgumentException("h cannot be zero");
		}
		double con = 1.4;
		double con2 = con * con;
		int nTab = 10;
		double safe = 2.0;
		double error = Double.MAX_VALUE;
		double result = 0.0, errorTotal, factor;
		double[][] m = new double[nTab][nTab];

		m[0][0] = (fn.evaluateAt(x + h) - fn.evaluateAt(x - h)) / (2.0 * h);
		for(int i = 1; i < nTab; ++i) {
			h /= con;
			m[0][i] = (fn.evaluateAt(x + h) - fn.evaluateAt(x - h)) / (2.0 * h);
			factor = con2;
			for(int j = 1; j <= i; ++j) {
				m[j][i] = (m[j - 1][i] * factor - m[j - 1][i - 1]) / (factor - 1.0);
				factor *= con2;
				errorTotal = Math.max(Math.abs(m[j][i] - m[j - 1][i]), Math.abs(m[j][i] - m[j - 1][i - 1]));
				if(errorTotal <= error) {
					error = errorTotal;
					result = m[j][i];
				}
			}
			if(Math.abs(m[i][i] - m[i - 1][i - 1]) >= safe * error) {
				break;
			}
		}
		return result;
	}

	public static void main(String[] args) {
		System.out.println(ridders(Math::tan, 1, 0.1));
	}
}
