package com.wildbitsfoundry.etk4j.math.interpolation;

import java.util.Arrays;

public class CardinalSpline extends Spline {

	private CardinalSpline(double[] x, double[] y, double[] dydx) {
		super(x, 3, y[0], y[y.length - 1]);
		final int n = _x.length;
		_coefs = new double[(n - 1) * 4]; // 4 coefficients and n - 1 segments
		for (int i = 0, j = 0, k = 0; i < n - 1; ++i, ++j, ++k) {
			double hx = _x[i + 1] - _x[i];
			double m0 = dydx[k];
			double m1 = dydx[++k];
			double a = (2 * y[i] + m0 - 2 * y[i + 1] + m1) / (hx * hx * hx);
			double b = (-3 * y[i] - 2 * m0 + 3 * y[i + 1] - m1) / (hx * hx);
			double c = m0 / hx;
			double d = y[i];
			_coefs[j] = a;
			_coefs[++j] = b;
			_coefs[++j] = c;
			_coefs[++j] = d;
		}
	}

	public static CardinalSpline newCardinalSpline(double[] x, double[] y, double tau, double alpha) {
		return newCardinalSplineInPlace(Arrays.copyOf(x, x.length), y, tau, alpha);
	}
	
	public static CardinalSpline newCardinalSplineInPlace(double[] x, double[] y, double tau, double alpha) {
		Splines.checkXYDimensions(x, y);
		Splines.checkMinkXLength(x, 4);
		final int n = x.length;
		double[] d = new double[n + 2];
		double cp = 1 - tau;

		// left end point
		double p0 = 2.0 * y[0] - y[1];
		double dx = x[1] - x[0]; // assume same distance to projected point
		double dy = y[0] - p0;
		double dt0 = Math.pow(dx * dx + dy * dy, alpha);

		dx = x[1] - x[0];
		dy = y[1] - y[0];
		double dt1 = Math.pow(dx * dx + dy * dy, alpha);

		d[0] = cp * (y[0] - p0) / dt0 - (y[1] - p0) / (dt0 + dt1) + (y[1] - y[0]) / dt1;
		d[0] *= dt1;

		dx = x[2] - x[1];
		dy = y[2] - y[1];
		double dt2 = Math.pow(dx * dx + dy * dy, alpha);

		d[1] = cp * (y[1] - y[0]) / dt1 - (y[2] - y[0]) / (dt1 + dt2) + (y[2] - y[1]) / dt2;
		d[1] *= dt1;

		// internal knots
		for (int i = 1; i < n - 1; ++i) {
			dx = x[i] - x[i - 1];
			dy = y[i] - y[i - 1];
			dt0 = Math.pow(dx * dx + dy * dy, alpha);

			dx = x[i + 1] - x[i];
			dy = y[i + 1] - y[i];
			dt1 = Math.pow(dx * dx + dy * dy, alpha);
			d[i + 1] = cp * (y[i] - y[i - 1]) / dt0 - (y[i + 1] - y[i - 1]) / (dt0 + dt1) + (y[i + 1] - y[i]) / dt1;
			d[i + 1] *= dt1;
			++i;

			dx = x[i + 1] - x[i];
			dy = y[i + 1] - y[i];
			dt2 = Math.pow(dx * dx + dy * dy, alpha);
			d[i + 1] = cp * (y[i] - y[i - 1]) / dt1 - (y[i + 1] - y[i - 1]) / (dt1 + dt2) + (y[i + 1] - y[i]) / dt2;
			d[i + 1] *= dt1;
		}

		// right end point
		final int dl = n + 2;
		dx = x[n - 2] - x[n - 3];
		dy = y[n - 3] - y[n - 3];
		dt0 = Math.pow(dx * dx + dy * dy, alpha);

		dx = x[n - 1] - x[n - 2];
		dy = y[n - 1] - y[n - 2];
		dt1 = Math.pow(dx * dx + dy * dy, alpha);

		d[dl - 2] = cp * (y[n - 2] - y[n - 3]) / dt0 - (y[n - 1] - y[n - 3]) / (dt0 + dt1)
				+ (y[n - 1] - y[n - 2]) / dt1;
		d[dl - 2] *= dt1;

		double p3 = 2.0 * y[n - 1] - y[n - 2];
		dx = x[n - 1] - x[n - 2]; // assume same distance to projected point
		dy = p3 - y[n - 1];
		dt2 = Math.pow(dx * dx + dy * dy, alpha);

		d[dl - 1] = cp * (y[n - 1] - y[n - 2]) / dt1 - (p3 - y[n - 2]) / (dt1 + dt2) + (p3 - y[n - 1]) / dt2;
		d[dl - 1] *= dt1;

		return new CardinalSpline(x, y, d);
	}
	

	public static CardinalSpline newCatmullRomSpline(double[] x, double[] y) {
		return newUniformCatmullRomSpline(x, y);
	}

	public static CardinalSpline newCentripetalCatmullRomSpline(double[] x, double[] y) {
		return CardinalSpline.newCardinalSpline(x, y, 0.0, 0.5);
	}

	public static CardinalSpline newUniformCatmullRomSpline(double[] x, double[] y) {
		return CardinalSpline.newCardinalSpline(x, y, 0.0, 0.0);
	}

	public static CardinalSpline newChordalCatmullRomSpline(double[] x, double[] y) {
		return CardinalSpline.newCardinalSpline(x, y, 0.0, 1.0);
	}

	// @Override
	// public String toString() {
	// StringBuilder sb = new StringBuilder();
	// for (int i = 0; i < _x.length - 1; ++i) {
	// double m0 = _dydx[i];
	// double m1 = _dydx[i + 1];
	// double a = 2 * _y[i] + m0 - 2 * _y[i + 1] + m1;
	// double b = -3 * _y[i] - 2 * m0 + 3 * _y[i + 1] - m1;
	// double c = m0;
	// double d = _y[i];
	//
	// sb.append("S").append(i + 1).append("(x) = ")
	// .append(a != 0d ? String.format("%.4f * (x - %.4f)^3", a, _x[i]) : "")
	// .append(b != 0d ? String.format(" + %.4f * (x - %.4f)^2", b, _x[i]) : "")
	// .append(c != 0d ? String.format(" + %.4f * (x - %.4f)", c, _x[i]) : "")
	// .append(d != 0d ? String.format(" + %.4f", d, _x[i]) : "")
	// .append(System.lineSeparator());
	// }
	// sb.setLength(Math.max(sb.length() - System.lineSeparator().length(), 0));
	// return sb.toString().replace("+ -", "- ");
	// }


	@Override
	public double differentiate(double x) {
		throw new RuntimeException("Method not implemented yet");
	}

	private final double a = 1.0 / 2.0, b = 1.0 / 3.0, c = 1.0 / 4.0;

	@Override
	protected double evaluateAntiDerivativeAt(int i, double t) {
		i <<= 2;
		return t * (_coefs[i + 3] + t * (_coefs[i + 2] * a + t * (_coefs[i + 1] * b + t * _coefs[i] * c)));
	}
}
