package com.wildbitsfoundry.etk4j.math.interpolation;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.polynomials.Polynomials;
import com.wildbitsfoundry.etk4j.util.ArrayUtils;

public class CubicSpline extends Spline {

	private static class TridiagonalLDLTSystem {
		public double[] D;
		public double[] SD;
		public double[] dydx;

		public TridiagonalLDLTSystem(int n) {
			SD = new double[n - 1];
			D = new double[n];
			dydx = new double[n];
		}

		public double[] solve() {
			tridiagonalLDLTSolve(D, SD, dydx);
			return dydx;
		}

		public double[] solve(int length, int offset) {
			tridiagonalLDLTSolve(D, SD, dydx, length, offset);
			return dydx;
		}

		private static void tridiagonalLDLTSolve(double[] D, double[] L, double[] b, int length, int offset) {
			int n = length;
			for (int i = offset; i < n - 1; ++i) {
				double ui = L[i];
				L[i] /= D[i];
				D[i + 1] -= ui * L[i];
				b[i + 1] -= L[i] * b[i];
			}

			b[n - 1] /= D[n - 1];
			for (int i = n - 2; i >= offset; --i) {
				b[i] = b[i] / D[i] - L[i] * b[i + 1];
			}
		}

		private static void tridiagonalLDLTSolve(double[] D, double[] L, double[] b) {
			tridiagonalLDLTSolve(D, L, b, b.length, 0);
		}
	}

	private interface CubicSplineBuilder {
		double[] build(TridiagonalLDLTSystem T, double[] x, double[] y, int n);
	}

	private CubicSpline(double[] x, double[] y, double[] dydx) {
		super(x, 4, y[0], y[y.length - 1]);

		final int n = _x.length;
		// compute coefficients
		_coefs = new double[(n - 1) * 4]; // 4 coefficients and n - 1 segments
		for (int i = 0, j = 0; i < n - 1; ++i, ++j) {
			double hx = _x[i + 1] - _x[i];
			if (hx <= 0.0) {
				throw new IllegalArgumentException("x must be monotonically increasing");
			}
			double m0 = dydx[i] * hx;
			double m1 = dydx[i + 1] * hx;
			double a = (2 * y[i] + m0 - 2 * y[i + 1] + m1) / (hx * hx * hx);
			double b = (-3 * y[i] - 2 * m0 + 3 * y[i + 1] - m1) / (hx * hx);
			double c = dydx[i];
			double d = y[i];
			_coefs[j] = a;
			_coefs[++j] = b;
			_coefs[++j] = c;
			_coefs[++j] = d;
		}
		double d0 = dydx[0];
		double dn = dydx[dydx.length - 1];
		super.setEndSlopes(d0, dn);
	}

	private CubicSpline(double[] x, double[] y, double[] coefs, double d0, double dn) {
		super(x, 4, y[0], y[y.length - 1]);
		_coefs = coefs;
		super.setEndSlopes(d0, dn);
	}

	@Override
	public void setExtrapolationMethod(ExtrapolationMethod method) {
		super.setExtrapolationMethod(method);
	}

	public static CubicSpline newCubicSpline(double[] x, double[] y) {
		return newNotAKnotSpline(x, y);
	}

	public static CubicSpline newCubicSplineInPlace(double[] x, double[] y) {
		return newNotAKnotSplineInPlace(x, y);
	}

	public static CubicSpline newNaturalSpline(double[] x, double[] y) {
		return newNaturalSplineInPlace(Arrays.copyOf(x, x.length), y);
	}

	private static CubicSpline buildSpline(double[] x, double[] y, int minLength, CubicSplineBuilder builder) {
		Splines.checkXYDimensions(x, y);
		Splines.checkMinkXLength(x, minLength);
		Polynomials.polyfit(x, y, 2);
		// double[] coefficients = parabola.GetCoefficients();
		// double a = coefficients[0];
		// double b = coefficients[1] - 2 * a * x[0];
		// double c = coefficients[2] - b * x[0] + a * x[0] * x[0];
		final int n = x.length;
		TridiagonalLDLTSystem T = setupC2Spline(x, y);
		return new CubicSpline(x, y, builder.build(T, x, y, n));
	}

	public static double[] parabolaFit(double x0, double x1, double x2, double y0, double y1, double y2) {
		double denom = (x0 - x1) * (x0 - x2) * (x1 - x2);
		double a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / denom;
		double b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / denom;
		double c = (x1 * x2 * (x1 - x2) * y0 + x2 * x0 * (x2 - x0) * y1 + x0 * x1 * (x0 - x1) * y2) / denom;

		return new double[] { a, b, c };

	}

	public static CubicSpline newNaturalSplineInPlace(double[] x, double[] y) {
		CubicSplineBuilder builder = new CubicSplineBuilder() {
			@Override
			public double[] build(TridiagonalLDLTSystem T, double[] x, double[] y, int n) {
				T.D[0] = 2.0 * T.SD[0];
				double r0 = (y[1] - y[0]) * T.SD[0] * T.SD[0];
				T.dydx[0] = 3.0 * r0;
				T.D[n - 1] = 2.0 * T.SD[n - 2];
				r0 = (y[n - 1] - y[n - 2]) * T.SD[n - 2] * T.SD[n - 2];
				T.dydx[n - 1] = 3.0 * r0;
				return T.solve();
			}
		};
		return buildSpline(x, y, 2, builder);
	}

	public static CubicSpline newParabolicallyTerminatedSpline(double[] x, double[] y) {

		return newParabolicallyTerminatedSplineInPlace(Arrays.copyOf(x, x.length), y);
	}

	public static CubicSpline newParabolicallyTerminatedSplineInPlace(double[] x, double[] y) {
		CubicSplineBuilder builder = new CubicSplineBuilder() {
			@Override
			public double[] build(TridiagonalLDLTSystem T, double[] x, double[] y, int n) {
				T.D[0] = T.SD[0];
				double r0 = (y[1] - y[0]) * T.SD[0] * T.SD[0];
				T.dydx[0] = 2.0 * r0;
				T.D[n - 1] = T.SD[n - 2];
				r0 = (y[n - 1] - y[n - 2]) * T.SD[n - 2] * T.SD[n - 2];
				T.dydx[n - 1] = 2.0 * r0;
				return T.solve();
			}
		};
		return buildSpline(x, y, 2, builder);
	}

	public static CubicSpline newClampedSpline(double[] x, double[] y, double d0, double dn) {
		return newClampedSplineInPlace(Arrays.copyOf(x, x.length), y, d0, dn);
	}

	public static CubicSpline newClampedSplineInPlace(double[] x, double[] y, double d0, double dn) {
		CubicSplineBuilder builder = new CubicSplineBuilder() {
			@Override
			public double[] build(TridiagonalLDLTSystem T, double[] x, double[] y, int n) {
				T.dydx[0] = d0;
				T.dydx[T.dydx.length - 1] = dn;
				T.dydx[1] = T.dydx[1] - T.dydx[0] * T.SD[0];
				T.dydx[n - 2] = T.dydx[n - 2] - T.dydx[n - 1] * T.SD[n - 2];
				return T.solve(T.dydx.length - 1, 1);
			}
		};
		return buildSpline(x, y, 2, builder);
	}

	public static CubicSpline newNotAKnotSpline(double[] x, double[] y) {
		return newNotAKnotSplineInPlace(Arrays.copyOf(x, x.length), y);
	}

	public static CubicSpline newNotAKnotSplineInPlace(double[] x, double[] y) {
		CubicSplineBuilder builder = new CubicSplineBuilder() {
			@Override
			public double[] build(TridiagonalLDLTSystem T, double[] x, double[] y, int n) {
				double a = T.SD[1] / T.SD[0];
				double b = 1 / (1 + a);
				T.D[0] = T.SD[0] * b;
				double r0 = (y[1] - y[0]) * T.SD[0] * T.SD[0];
				double r1 = (y[2] - y[1]) * T.SD[1] * T.SD[1];
				T.dydx[0] = ((3.0 * a + 2.0) * r0 + a * r1) * b * b;

				a = T.SD[n - 3] / T.SD[n - 2];
				b = 1 / (1 + a);
				T.D[n - 1] = T.SD[n - 2] * b;
				r0 = (y[n - 1] - y[n - 2]) * T.SD[n - 2] * T.SD[n - 2];

				r1 = (y[n - 2] - y[n - 3]) * T.SD[n - 3] * T.SD[n - 3];
				T.dydx[n - 1] = ((3.0 * a + 2.0) * r0 + a * r1) * b * b;
				return T.solve();
			}
		};
		return buildSpline(x, y, 4, builder);
	}

	public static CubicSpline newAkimaSpline(double[] x, double[] y) {
		return newAkimaSpline(x, y, 1e-12);
	}

	public static CubicSpline newAkimaSplineInPlace(double[] x, double[] y) {
		return newAkimaSplineInPlace(x, y, 1e-12);
	}

	public static CubicSpline newAkimaSpline(double[] x, double[] y, double ep) {
		return newAkimaSplineInPlace(Arrays.copyOf(x, x.length), y, ep);
	}

	public static CubicSpline newAkimaSplineInPlace(double[] x, double[] y, double ep) {
		Splines.checkXYDimensions(x, y);
		Splines.checkMinkXLength(x, 5);
		final int n = x.length;
		double[] d = new double[n];
		double[] t = new double[n + 3];

		// inner knots
		for (int i = 0; i < n - 1; ++i) {
			t[i + 2] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
		}
		// end points
		t[n + 1] = 2.0 * t[n] - t[n - 1];
		t[n + 2] = 2.0 * t[n + 1] - t[n];
		t[1] = 2.0 * t[2] - t[3];
		t[0] = 2.0 * t[1] - t[2];

		for (int i = 0; i < n; ++i) {
			double c0 = Math.abs(t[i + 3] - t[i + 2]);
			double c1 = Math.abs(t[i + 1] - t[i]);
			double c2 = c0 + c1;
			if (c2 > ep) {
				d[i] = (c0 * t[i + 1] + c1 * t[i + 2]) / c2;
			} else {
				d[i] = 0.5 * (t[i + 2] + t[i + 1]);
			}
		}
		return new CubicSpline(x, y, d);
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0, j = 0; i < _x.length - 1; ++i, ++j) {
			double a = _coefs[j];
			double b = _coefs[++j];
			double c = _coefs[++j];
			double d = _coefs[++j];

			sb.append("S").append(i + 1).append("(x) = ")
					.append(a != 0d ? String.format("%.4f * (x - %.4f)^3", a, _x[i]) : "")
					.append(b != 0d ? String.format(" + %.4f * (x - %.4f)^2", b, _x[i]) : "")
					.append(c != 0d ? String.format(" + %.4f * (x - %.4f)", c, _x[i]) : "")
					.append(d != 0d ? String.format(" + %.4f", d, _x[i]) : "").append(System.lineSeparator());
		}
		sb.setLength(Math.max(sb.length() - System.lineSeparator().length(), 0));
		return sb.toString().replace("+ -", "- ").replace("- -", "+ ");
	}

	// move to examples
	public static void main(String[] args) {
		double[] x = new double[] { 1, 2, 3, 4 };
		double[] y = new double[] { 5, 6, 7, 8 };
		CubicSpline cs = newNotAKnotSpline(x, y);

		double[] xi = ArrayUtils.linsteps(1.0, 4.0, 0.1);
		for (double xii : xi) {
			System.out.printf("y(%.4f) = %.4f%n", xii, cs.evaluateAt(xii));
		}

	}

	private static TridiagonalLDLTSystem setupC2Spline(double[] x, double[] y) {
		final int n = x.length - 1;
		TridiagonalLDLTSystem T = new TridiagonalLDLTSystem(n + 1);
		double[] D = T.D;
		double[] SD = T.SD;
		double[] dydx = T.dydx;

		SD[0] = 1.0 / (x[1] - x[0]);
		double r1 = (y[1] - y[0]) * SD[0] * SD[0];
		for (int i = 1; i < n; ++i) {
			SD[i] = 1.0 / (x[i + 1] - x[i]);
			double r0 = r1;

			r1 = (y[i + 1] - y[i]) * SD[i] * SD[i];
			D[i] = 2.0 * (SD[i - 1] + SD[i]);
			dydx[i] = 3.0 * (r0 + r1);
		}
		return T;
	}

	@Override
	protected double getValueAt(int i, double x) {
		double t = x - _x[i];
		i <<= 2;
		return _coefs[i + 3] + t * (_coefs[i + 2] + t * (_coefs[i + 1] + t * _coefs[i]));
	}

	@Override
	public double evaluateDerivativeAt(int i, double t) {
		i <<= 2;
		return _coefs[i + 2] + t * (_coefs[i + 1] + t * _coefs[i]);
	}

	private final double a = 1.0 / 2.0, b = 1.0 / 3.0, c = 1.0 / 4.0;

	@Override
	protected double evaluateAntiDerivativeAt(int i, double t) {
		i <<= 2;
		return t * (_coefs[i + 3] + t * (_coefs[i + 2] * a + t * (_coefs[i + 1] * b + t * _coefs[i] * c)));
	}
}