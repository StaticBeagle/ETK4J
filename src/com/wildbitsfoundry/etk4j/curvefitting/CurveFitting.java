package com.wildbitsfoundry.etk4j.curvefitting;

import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;

public final class CurveFitting {
	private CurveFitting() {
	}

	public static double[] line(double x0, double x1, double y0, double y1) {
		double m = (y1 - y0) / (x1 - x0);
		return new double[] { m, y0 - m * x0 };
	}

	public static double[] parabola(double x0, double x1, double x2, double y0, double y1, double y2) {
		double denom = (x0 - x1) * (x0 - x2) * (x1 - x2);
		double a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / denom;
		double b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / denom;
		double c = (x1 * x2 * (x1 - x2) * y0 + x2 * x0 * (x2 - x0) * y1 + x0 * x1 * (x0 - x1) * y2) / denom;

		return new double[] { a, b, c };
	}

	public double[] polynomial(double[] x, double[] y, int n) {
		return Polynomial.polyfit(x, y, n).getCoefficients();
	}

	public double[] exponential(double[] x, double[] y) {
		double a = 0.0, b = 0.0, c = 0.0, d = 0.0, r0 = 0.0, r1 = 0.0;
		for (int i = 0; i < x.length; i++) {
			a += y[i];
			d = x[i] * y[i];
			b += d;
			c += d * x[i];
			d = y[i] * Math.log(y[i]);
			r0 += d;
			r1 += x[i] * d;
		}
		double det = 1 / (a * c - b * b);
		double alpha = Math.exp((c * r0 - b * r1) * det);
		double beta = (a * r1 - b * r0) * det;
		return new double[] { alpha, beta };
	}

}
