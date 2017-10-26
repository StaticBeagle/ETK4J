package com.wildbitsfoundry.etk4j.math.interpolation;

public class QuadraticSpline extends Spline {

	private static class TridiagonalSystem {
		public double[] L; // sub-diagonal
		public double[] D; // diagonal
		public double[] U; // super-diagonal
		public double[] dydx;

		public TridiagonalSystem(int n) {
			L = new double[n - 1];
			D = new double[n];
			U = new double[n - 1];
			dydx = new double[n];
		}

		public double[] solve() {
			int n = dydx.length;
			for (int i = 0; i < n - 1; ++i) {
				L[i] /= D[i];
				D[i + 1] -= U[i] * L[i];
				dydx[i + 1] -= L[i] * dydx[i];
			}

			dydx[n - 1] /= D[n - 1];
			for (int i = n - 2; i >= 0; --i) {
				dydx[i] = (dydx[i] - U[i] * dydx[i + 1]) / D[i];
			}
			return dydx;
		}
	}

	protected QuadraticSpline(double[] x, int order, double y0, double yn) {
		super(x, 3, y0, yn);
	}

	public static QuadraticSpline newQuadraticSpline(double[] x, double[] y) {


		return null;
	}

	public static void main(String[] args) {
		TridiagonalSystem T = new TridiagonalSystem(3);
		T.D = new double[] { 1, 1, 1 };
		T.L = new double[] { 1, 1 };
		T.U = new double[] { 0, 0 };

		T.dydx = new double[] { 0, 2, 4 };

		T.solve();
	}

	@Override
	protected double evaluateDerivativeAt(int i, double t) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	protected double evaluateAntiDerivativeAt(int i, double t) {
		// TODO Auto-generated method stub
		return 0;
	}

}
