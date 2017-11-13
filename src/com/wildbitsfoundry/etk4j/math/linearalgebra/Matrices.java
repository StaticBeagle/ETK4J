package com.wildbitsfoundry.etk4j.math.linearalgebra;

public final class Matrices {
	private Matrices() {
	}
	
	public static Matrix fwdSubsSolve(Matrix L, Matrix B) {
		final int nx = B.getColumnCount();
		final int m = L.getRowCount();
		final int n = L.getColumnCount();
		double[] t = L.getArray();
		double[] X = B.getArrayCopy();

		for(int j = 0; j < nx; ++j) {
			for(int i = 0; i < m; ++i) {
				for(int k = 0; k < i; ++k) {
					X[i * nx + j] -= X[k * nx + j] * t[i * n + k];
				}
				X[i * nx + j] /= t[i * n + i];
			}			
		}
		return new Matrix(X, m, nx);
	}

//	public static double[] backSubtitution(Matrix U, double[] b) {
//		final int length = b.length;
//		double[] x = new double[length];
//		for (int i = length - 1, index; i >= 0; --i) {
//			x[i] = b[i];
//			index = i * U._cols;
//			for (int j = i + 1; j < length; ++j) {
//				x[i] -= U._data[index + j] * x[j];
//			}
//			x[i] /= U._data[index + i];
//		}
//		return x;
//	}
	
	public static void tridiagonalLDLTSolve(double[] D, double[] L, double[] b, int length, int offset) {
		int n = length;
		for(int i = offset; i < n - 1; ++i) {
			double ui = L[i];
			L[i] /= D[i];
			D[i + 1] -= ui * L[i];
			b[i + 1] -= L[i] * b[i];
		}

		b[n - 1] /= D[n - 1];
		for(int i = n - 2; i >= offset; --i) {
			b[i] = b[i] / D[i] - L[i] * b[i + 1];
		}
	}
	
	public static void tridiagonalLDLTSolve(double[] D, double[] L, double[] b) {
		tridiagonalLDLTSolve(D, L, b, b.length, 0);
	}

	public static Matrix Vandermonde(double[] x, int rows, int cols) {
		double[][] V = new double[rows][cols];
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				V[i][j] = Math.pow(x[i], j);
			}
		}
		return new Matrix(V);
	}

	public static Matrix Identity(int rows, int cols) {
		double[] data = new double[rows * cols];
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				if (i == j) {
					data[i * cols + j] = 1.0;
				}
			}
		}
		return new Matrix(data, rows, cols);
	}

	public static Matrix Identity(int n) {
		return Matrices.Identity(n, n);
	}

	public static Matrix Companion(double[] coefs, int n) {
		// Construct the companion matrix
		Matrix c = new Matrix(n, n);

		double a = 1.0 / coefs[0];
		for (int i = 0; i < n; i++) {
			c.set(0, n - 1 - i, -coefs[n - i] * a);
		}
		for (int i = 1; i < n; i++) {
			c.set(i, i - 1, 1);
		}
		return c;
	}
}
