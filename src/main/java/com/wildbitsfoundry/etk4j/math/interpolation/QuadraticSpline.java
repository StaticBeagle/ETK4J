package com.wildbitsfoundry.etk4j.math.interpolation;

import java.util.Arrays;

public class QuadraticSpline extends Spline {
	
	private static final double P5 = 0.5, P33 = 1.0 / 3.0;

	private static class TridiagonalSystem {
		public double[] L; // Sub-diagonal
		public double[] D; // Diagonal
		public double[] U; // Super-diagonal
		public double[] b; // Solution vector

		public TridiagonalSystem(int n) {
			L = new double[n - 1];
			D = new double[n];
			U = new double[n - 1];
			b = new double[n];
		}

		public double[] solve() {
			int n = b.length;
			for (int i = 0; i < n - 1; ++i) {
				L[i] /= D[i];
				D[i + 1] -= U[i] * L[i];
				b[i + 1] -= L[i] * b[i];
			}

			b[n - 1] /= D[n - 1];
			for (int i = n - 2; i >= 0; --i) {
				b[i] = (b[i] - U[i] * b[i + 1]) / D[i];
			}
			return b;
		}
	}
	
	protected QuadraticSpline(double[] x, double[] y, double[] dydx) {
		super(x, 3);

		final int n = _x.length;
		// compute coefficients
		_coefs = new double[(n - 1) * 3]; // 3 coefficients and n - 1 segments
		for (int i = 0, j = 0; i < n - 1; ++i, ++j) {
			double hx = _x[i + 1] - _x[i];
			if (hx <= 0.0) {
				throw new IllegalArgumentException("x must be monotonically increasing");
			}
			double a = 0.5 * (dydx[i + 1] - dydx[i]) / hx;
			double b = dydx[i];
			double c = y[i];
			_coefs[j] = a;
			_coefs[++j] = b;
			_coefs[++j] = c;
		}
	}

	public static QuadraticSpline newNaturalSpline(double[] x, double[] y) {
		return newNaturalSplineInPlace(Arrays.copyOf(x, x.length), y);
	}

	public static QuadraticSpline newNaturalSplineInPlace(double[] x, double[] y) {

		TridiagonalSystem T = setupSpline(x, y);
		// Natural conditions
		T.b[0] = 0.0;
		T.D[0] = 1.0;

		return new QuadraticSpline(x, y, T.solve());
	}

	public static QuadraticSpline newClampedSpline(double[] x, double[] y, double d0) {
		return newClampedSplineInPlace(Arrays.copyOf(x, x.length), y, d0);
	}

	public static QuadraticSpline newClampedSplineInPlace(double[] x, double[] y, double d0) {

		TridiagonalSystem T = setupSpline(x, y);
		// Natural conditions
		T.b[0] = d0;
		T.D[0] = 1.0;

		return new QuadraticSpline(x, y, T.solve());
	}

	private static TridiagonalSystem setupSpline(double[] x, double[] y) {
		final int n = x.length;
		TridiagonalSystem T = new TridiagonalSystem(n);

		// U is always zero
		// L and D are always one
		Arrays.fill(T.L, 1.0);
		Arrays.fill(T.D, 1.0);

		for (int i = 1; i < n; ++i) {
			double hx = x[i] - x[i - 1];
			T.b[i] = 2 * (y[i] - y[i - 1]) / hx;
		}
		return T;
	}

	

	@Override
	public double evaluateSegmentAt(int i, double x) {
		double t = x - _x[i];
		i *= 3;
		return _coefs[i + 2] + t * (_coefs[i + 1] + t * _coefs[i]);
	}

	@Override
	public double evaluateDerivativeAt(int i, double t) {
		i *= 3;
		return _coefs[i + 1] + t * 2.0 * _coefs[i];
	}

	@Override
	public double evaluateAntiDerivativeAt(int i, double t) {
		i *= 3;
		return t * (_coefs[i + 2] + t * (_coefs[i + 1] * P5 + t * _coefs[i] * P33));
	}

	public static void main(String[] args) {
		double[] x = { -1, 0, 1 };
		double[] y = { 0, 1, 3 };
		QuadraticSpline qs = newNaturalSpline(x, y);
		
		System.out.println(qs.evaluateAt(-0.5));
		System.out.println(qs.evaluateAt(0.5));
		
		x = new double[] {0, 1, 2, 3, 4, 5, 6, 7, 8 };
		y = new double[] {0, 1, 4, 9, 16, 25, 36, 49, 64 };
		
		QuadraticSpline cs = QuadraticSpline.newNaturalSpline(x, y);
		System.out.println(cs.differentiate(2.0));
		System.out.println(cs.evaluateAntiDerivativeAt(0, 2.0));
		System.out.println(cs.integrate(3.0));
		System.out.println(cs.integrate(1, 3.0));
	}
}
