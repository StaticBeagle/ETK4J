package com.wildbitsfoundry.etk4j.math.interpolation;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkMinXLength;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkXYDimensions;

public class CubicSpline extends Spline {
	
	private static final double P5 = 0.5, P33 = 1.0 / 3.0, P25 = 0.25;

	private static class TridiagonalLDLTSystem {
		public double[] D;
		public double[] L;
		public double[] b;

		public TridiagonalLDLTSystem(int n) {
			L = new double[n - 1];
			D = new double[n];
			b = new double[n];
		}

		public double[] solve() {
			return this.solve(b.length, 0);
		}

		public double[] solve(int length, int offset) {
			for (int i = offset; i < length - 1; ++i) {
				double ui = L[i];
				L[i] /= D[i];
				D[i + 1] -= ui * L[i];
				b[i + 1] -= L[i] * b[i];
			}

			b[length - 1] /= D[length - 1];
			for (int i = length - 2; i >= offset; --i) {
				b[i] = b[i] / D[i] - L[i] * b[i + 1];
			}
			return b;
		}
	}

	private CubicSpline(double[] x, double[] y, double[] dydx) {
		super(x, 4);
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
	}

	private CubicSpline(double[] x, double[] y, double[] coefs, double d0, double dn) {
		super(x, 4);
		_coefs = coefs;
	}

	private static TridiagonalLDLTSystem setupLDLT(double[] x, double[] y) {
		checkXYDimensions(x, y);
		final int n = x.length;
		TridiagonalLDLTSystem T = new TridiagonalLDLTSystem(n);
		double[] D = T.D;
		double[] SD = T.L;
		double[] dydx = T.b;

		SD[0] = 1.0 / (x[1] - x[0]);
		double r1 = (y[1] - y[0]) * SD[0] * SD[0];
		for (int i = 1; i < n - 1; ++i) {
			SD[i] = 1.0 / (x[i + 1] - x[i]);
			double r0 = r1;

			r1 = (y[i + 1] - y[i]) * SD[i] * SD[i];
			D[i] = 2.0 * (SD[i - 1] + SD[i]);
			dydx[i] = 3.0 * (r0 + r1);
		}
		return T;
	}

	public static CubicSpline newHermiteSplineInPlace(double[] x, double[] y, double[] dydx) {
		checkXYDimensions(x, y);
		checkXYDimensions(x, dydx);
		return new CubicSpline(x, y, dydx);
	}

	public static CubicSpline newHermiteSpline(double[] x, double[] y, double[] dydx) {
		return newHermiteSplineInPlace(Arrays.copyOf(x, x.length), y, dydx);
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

	public static double[] computeNaturalSplineDerivatives(double[] x, double[] y) {
		checkMinXLength(x, 2);
		TridiagonalLDLTSystem T = setupLDLT(x, y);
		final int n = x.length;
		T.D[0] = 2.0 * T.L[0];
		double r0 = (y[1] - y[0]) * T.L[0] * T.L[0];
		T.b[0] = 3.0 * r0;
		T.D[n - 1] = 2.0 * T.L[n - 2];
		r0 = (y[n - 1] - y[n - 2]) * T.L[n - 2] * T.L[n - 2];
		T.b[n - 1] = 3.0 * r0;
		return T.solve();
	}

	public static CubicSpline newNaturalSplineInPlace(double[] x, double[] y) {
		double[] dydx = computeNaturalSplineDerivatives(x, y);
		return new CubicSpline(x, y, dydx);
	}

	public static CubicSpline newParabolicallyTerminatedSpline(double[] x, double[] y) {
		return newParabolicallyTerminatedSplineInPlace(Arrays.copyOf(x, x.length), y);
	}

	public static double[] computeParabolicallyTerminatedSplineDerivatives(double[] x, double[] y) {
		checkMinXLength(x, 2);
		TridiagonalLDLTSystem T = setupLDLT(x, y);
		final int n = x.length;
		T.D[0] = T.L[0];
		double r0 = (y[1] - y[0]) * T.L[0] * T.L[0];
		T.b[0] = 2.0 * r0;
		T.D[n - 1] = T.L[n - 2];
		r0 = (y[n - 1] - y[n - 2]) * T.L[n - 2] * T.L[n - 2];
		T.b[n - 1] = 2.0 * r0;
		return T.solve();
	}

	public static CubicSpline newParabolicallyTerminatedSplineInPlace(double[] x, double[] y) {
		double[] dydx = computeParabolicallyTerminatedSplineDerivatives(x, y);
		return new CubicSpline(x, y, dydx);
	}

	public static double[] computeClampedSplineDerivatives(double[] x, double[] y, double d0, double dn) {
		checkMinXLength(x, 2);
		TridiagonalLDLTSystem T = setupLDLT(x, y);
		final int n = x.length;
		T.b[0] = d0;
		T.b[T.b.length - 1] = dn;
		T.b[1] = T.b[1] - T.b[0] * T.L[0];
		T.b[n - 2] = T.b[n - 2] - T.b[n - 1] * T.L[n - 2];
		return T.solve(T.b.length - 1, 1);
	}

	public static CubicSpline newClampedSpline(double[] x, double[] y, double d0, double dn) {
		return newClampedSplineInPlace(Arrays.copyOf(x, x.length), y, d0, dn);
	}

	public static CubicSpline newClampedSplineInPlace(double[] x, double[] y, double d0, double dn) {
		double[] dydx = computeClampedSplineDerivatives(x, y, d0, dn);
		return new CubicSpline(x, y, dydx);
	}

	public static CubicSpline newNotAKnotSpline(double[] x, double[] y) {
		return newNotAKnotSplineInPlace(Arrays.copyOf(x, x.length), y);
	}

	public static double[] computeNotAKnotSplineDerivatives(double[] x, double[] y) {
		checkMinXLength(x, 4);
		TridiagonalLDLTSystem T = setupLDLT(x, y);
		final int n = x.length;
		double a = T.L[1] / T.L[0];
		double b = 1 / (1 + a);
		T.D[0] = T.L[0] * b;
		double r0 = (y[1] - y[0]) * T.L[0] * T.L[0];
		double r1 = (y[2] - y[1]) * T.L[1] * T.L[1];
		T.b[0] = ((3.0 * a + 2.0) * r0 + a * r1) * b * b;

		a = T.L[n - 3] / T.L[n - 2];
		b = 1 / (1 + a);
		T.D[n - 1] = T.L[n - 2] * b;
		r0 = (y[n - 1] - y[n - 2]) * T.L[n - 2] * T.L[n - 2];

		r1 = (y[n - 2] - y[n - 3]) * T.L[n - 3] * T.L[n - 3];
		T.b[n - 1] = ((3.0 * a + 2.0) * r0 + a * r1) * b * b;
		return T.solve();
	}

	public static CubicSpline newNotAKnotSplineInPlace(double[] x, double[] y) {
		double[] dydx = computeNotAKnotSplineDerivatives(x, y);
		return new CubicSpline(x, y, dydx);
	}

	public static CubicSpline newAkimaSpline(double[] x, double[] y) {
		return newAkimaSpline(x, y, ConstantsETK.DOUBLE_EPS);
	}

	public static CubicSpline newAkimaSplineInPlace(double[] x, double[] y) {
		return newAkimaSplineInPlace(x, y, ConstantsETK.DOUBLE_EPS);
	}

	public static CubicSpline newAkimaSpline(double[] x, double[] y, double ep) {
		return newAkimaSplineInPlace(Arrays.copyOf(x, x.length), y, ep);
	}

	public static CubicSpline newAkimaSplineInPlace(double[] x, double[] y, double ep) {
		checkXYDimensions(x, y);
		checkMinXLength(x, 5);
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
	public double evaluateSegmentAt(int i, double x) {
		double t = x - _x[i];
		i <<= 2;
		return _coefs[i + 3] + t * (_coefs[i + 2] + t * (_coefs[i + 1] + t * _coefs[i]));
	}

	@Override
	public double evaluateDerivativeAt(int i, double t) {
		i <<= 2;
		return _coefs[i + 2] + t * (2 * _coefs[i + 1] + t * 3 * _coefs[i]);
	}

	@Override
	protected double evaluateAntiDerivativeAt(int i, double t) {
		i <<= 2;
		return t * (_coefs[i + 3] + t * (_coefs[i + 2] * P5 + t * (_coefs[i + 1] * P33 + t * _coefs[i] * P25)));
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
}