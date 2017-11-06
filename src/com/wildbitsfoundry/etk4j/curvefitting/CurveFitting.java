package com.wildbitsfoundry.etk4j.curvefitting;

import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.util.NumArrays;

public final class CurveFitting {
	private CurveFitting() {
	}

	public static double[] line(double x0, double x1, double y0, double y1) {
		double m = (y1 - y0) / (x1 - x0);
		return new double[] { m, y0 - m * x0 };
	}

	public static double[] linear(double[] x, double[] y) {
		if (x.length != y.length) {
			throw new IllegalArgumentException("x and y dimensions must match");
		}

		final int n = x.length;
		double xbar = NumArrays.mean(x);
		double ybar = NumArrays.mean(y);

		double sumx2 = NumArrays.norm(x);
		sumx2 *= sumx2;

		double ssxx = sumx2 - n * xbar * xbar;

		double xdoty = 0.0;
		for (int i = 0; i < n; ++i) {
			xdoty += x[i] * y[i];
		}
		double ssxy = xdoty - n * xbar * ybar;

		double b = ssxy / ssxx;
		double a = ybar - b * xbar;
		return new double[] { b, a };
	}

	public static double[] parabola(double x0, double x1, double x2, double y0, double y1, double y2) {
		double denom = (x0 - x1) * (x0 - x2) * (x1 - x2);
		double a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / denom;
		double b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / denom;
		double c = (x1 * x2 * (x1 - x2) * y0 + x2 * x0 * (x2 - x0) * y1 + x0 * x1 * (x0 - x1) * y2) / denom;

		return new double[] { a, b, c };
	}

	public double[] polynomial(double[] x, double[] y, int n) {
		return Polynomial.polyFit(x, y, n).getCoefficients();
	}

	/***
	 * Exponential least squares fit <br />
	 * Fits x and y to a function a * e^(b * x) in a least square sense
	 * 
	 * @param x
	 * @param y
	 * @return [a, b]
	 */
	public double[] exponential(double[] x, double[] y) {
		if (x.length != y.length) {
			throw new IllegalArgumentException("x and y dimensions must match");
		}
		final int n = x.length;
		double k0 = 0.0, k1 = 0.0, k2 = 0.0, k3 = 0.0, r0 = 0.0, r1 = 0.0;
		for (int i = 0; i < n; i++) {
			k0 += y[i];
			k3 = x[i] * y[i];
			k1 += k3;
			k2 += k3 * x[i];
			k3 = y[i] * Math.log(y[i]);
			r0 += k3;
			r1 += x[i] * k3;
		}
		double det = 1 / (k0 * k2 - k1 * k1);
		double a = Math.exp((k2 * r0 - k1 * r1) * det);
		double b = (k0 * r1 - k1 * r0) * det;
		return new double[] { a, b };
	}

	/***
	 * Logarithmic least squares fit Fits x and y to a function a + b * ln(x) in
	 * a least square sense
	 * 
	 * @param x
	 * @param y
	 * @return [a, b]
	 */
	public double[] logarithmic(double[] x, double[] y) {
		if (x.length != y.length) {
			throw new IllegalArgumentException("x and y dimensions must match");
		}
		final int n = x.length;
		double k0 = 0.0, k1 = 0.0, k2 = 0.0, k3 = 0.0, k4 = 0.0;
		for (int i = 0; i < n; i++) {
			k0 += y[i];
			k3 = Math.log(x[i]);
			k1 += k3;
			k2 += k3 * k3;
			k4 += y[i] * k3;
		}
		double b = (n * k4 - k0 * k1) / (n * k2 - k1 * k1);
		double a = (k0 - b * k1) / n;
		return new double[] { a, b };
	}
}
