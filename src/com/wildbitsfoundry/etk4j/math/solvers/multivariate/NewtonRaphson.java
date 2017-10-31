package com.wildbitsfoundry.etk4j.math.solvers.multivariate;

import com.wildbitsfoundry.etk4j.math.functions.MultivariateFunction;

public class NewtonRaphson {

	public static double[] Newton(double[] x0, MultivariateFunction[] functions, MultivariateFunction[][] jacobian) {
		int maxiter = 100;
		double tol = 1e-10;
		double maxval = 20e4;

		int rows = jacobian.length;
		int cols = jacobian[0].length;
		double[] xcurrent = x0;
		double[] xfinal = new double[cols];
		double[] functionValues = new double[cols];
		double[][] jacobianMatrixValues = new double[rows][cols];
		while (maxiter > 0) {
			evaluateJacobian(xcurrent, jacobianMatrixValues, jacobian);
			LU matrix = new LU(jacobianMatrixValues);
			if (Double.compare(Math.abs(matrix.det()), tol) < 0) {
				// bail out the matrix is singular
				return null; // <-- we'll think about this later
			}

			// evaluate functions at current x
			eval(xcurrent, functionValues, functions);
			double[] residuals = matrix.solve(functionValues);

			for (int i = 0; i < cols; i++) {
				xfinal[i] = xcurrent[i] - residuals[i];
			}
			eval(xfinal, functionValues, functions);

			boolean done = true;
			for (int i = 0; i < cols; i++) {
				done &= Double.compare(Math.abs(residuals[i]), tol) < 0;
			}
			if (done) {

				return xfinal;
			}

			boolean diverged = false;
			for (int i = 0; i < cols; i++) {
				diverged |= Double.compare(Math.abs(xcurrent[i]), maxval) > 0;
			}
			if (diverged) {
				return null;
			}

			maxiter--;
			xcurrent = xfinal;

		}
		return null;
	}

	public static class SolverResults {
		public boolean Converged;
		public double[] Results;
		public int NoIterations;
		public double Error;
		public String Status;
	}

	public static double[] Newton(final double[] x0, MultivariateFunction[] functions) {
		int maxiter = 100;
		double tol = 1e-15;
		double maxval = 20e4;
		double step = 0.001;
		double h = 1.0 / step;

		int rows = functions.length;
		int cols = x0.length;
		double[] xold = new double[cols];
		double[] xfinal = new double[cols];
		double[] residuals = new double[cols];
		double[] functionValues = new double[cols];
		double[] functionPrimeValues = new double[cols];
		double[][] jacobianMatrixValues = new double[rows][cols];

		System.arraycopy(x0, 0, xold, 0, rows);
		while (maxiter > 0) {
			// Compute Jacobian
			for (int i = 0; i < cols; i++) {
				xold[i] = xold[i] + step;
				eval(xold, functionPrimeValues, functions);
				xold[i] = xold[i] - step;
				eval(xold, functionValues, functions);
				for (int j = 0; j < rows; j++) {
					jacobianMatrixValues[j][i] = (functionPrimeValues[j] - functionValues[j]) * h;
				}
			}

			LU matrix = new LU(jacobianMatrixValues);
			if (Double.compare(Math.abs(matrix.det()), tol) < 0) {
				// bail out the matrix is singular
				return null; // <-- we'll think about this later
			}

			residuals = matrix.solve(functionValues);
			for (int i = 0; i < cols; i++) {
				xfinal[i] = xold[i] - residuals[i];
			}
			eval(xfinal, functionValues, functions);

			double relErrorNorm = 0.0;
			double errorNorm = 0.0;
			for (int i = 0; i < cols; i++) {
				// calculate the norm2 of the relative error
				relErrorNorm += (xfinal[i] - xold[i]) * (xfinal[i] - xold[i]);
				errorNorm += xfinal[i] * xfinal[i];
			}
			double error = Math.sqrt(relErrorNorm) / Math.sqrt(errorNorm);
			if (Double.compare(error, tol) <= 0) {
				return xfinal;
			}

			boolean diverged = false;
			for (int i = 0; i < cols; i++) {
				diverged |= Double.compare(Math.abs(xold[i]), maxval) > 0;
			}
			if (diverged) {
				return null;
			}

			maxiter--;
			System.arraycopy(xfinal, 0, xold, 0, rows);
		}
		return null;
	}

	public static void evaluateJacobian(double[] x, double[][] out, MultivariateFunction[][] jacobian) {
		int rows = out.length;
		int cols = out[0].length;

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				out[i][j] = jacobian[i][j].evaluateAt(x);
			}
		}
	}

	public static void computeJacobian(double[] x, double[][] out, MultivariateFunction[] functions, double step) {
		int rows = out.length;
		int cols = out[0].length;
		double h = 1.0 / step;
		double[] fxx = new double[rows];
		double[] fx = new double[rows];
		for (int i = 0; i < cols; i++) {
			x[i] = x[i] + step;
			eval(x, fxx, functions);
			x[i] = x[i] - step;
			eval(x, fx, functions);
			for (int j = 0; j < rows; j++) {
				out[j][i] = (fxx[j] - fx[j]) * h;
			}
		}
	}

	public static void eval(double[] x, double[] out, MultivariateFunction[] functions) {
		int length = x.length;
		for (int i = 0; i < length; i++) {
			out[i] = functions[i].evaluateAt(x);
		}
	}

	public static void main(String[] args) {
		MultivariateFunction f1 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return 3 * x[0] - Math.cos(x[1] * x[2]) - 1.5;
			}
		};

		MultivariateFunction f2 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return 4 * x[0] * x[0] - 625.0 * x[1] * x[1] + 2.0 * x[2] - 1.0;
			}
		};

		MultivariateFunction f3 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return 20 * x[2] + Math.exp(-x[0] * x[1]) + 9.0;
			}
		};

		// Jacobian
		MultivariateFunction df1x1 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return 3;
			}
		};

		MultivariateFunction df1x2 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return x[2] * Math.sin(x[1] * x[2]);
			}
		};

		MultivariateFunction df1x3 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return x[1] * Math.sin(x[1] * x[2]);
			}
		};

		MultivariateFunction df2x1 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return 8 * x[0];
			}
		};

		MultivariateFunction df2x2 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return -1250.0 * x[1];
			}
		};

		MultivariateFunction df2x3 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return 2.0;
			}
		};

		MultivariateFunction df3x1 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return -x[1] * Math.exp(x[0] * x[1]);
			}
		};

		MultivariateFunction df3x2 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return -x[0] * Math.exp(x[0] * x[1]);
			}
		};

		MultivariateFunction df3x3 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return 20.0;
			}
		};

		MultivariateFunction[][] Jacobian = new MultivariateFunction[][] { { df1x1, df1x2, df1x3 },
				{ df2x1, df2x2, df2x3 }, { df3x1, df3x2, df3x3 } };

		MultivariateFunction[] functions = new MultivariateFunction[] { f1, f2, f3 };
		double[] initialguess = new double[] { 1d, 2d, 3d };

		System.out.println("Solution using pre-computed Jacobian");
		for (double val :  Newton(initialguess, functions, Jacobian)) {
			System.out.println(val);
		}
		
		System.out.printf("%nSolution with auto-computed Jacobian%n");
		for (double val :  Newton(initialguess, functions)) {
			System.out.println(val);
		}
	}

	private static class LU {
		protected double[][] _data;
		protected int _cols;

		protected int _pivotsign;
		protected int[] _pivot;

		public LU(double[][] matrix) {
			int rows = matrix.length;
			int cols = _cols = matrix[0].length;
			_data = matrix;

			_pivot = new int[rows];
			for (int i = 0; i < rows; i++) {
				_pivot[i] = i;
			}
			_pivotsign = 1;

			double[] LUrow;
			double[] LUcol = new double[rows];

			// Begin the outer loop
			for (int j = 0; j < cols; j++) {
				// Copy the j-th column to localize references.
				for (int i = 0; i < rows; i++) {
					LUcol[i] = _data[i][j];
				}
				// Apply previous transformations
				for (int i = 0; i < rows; i++) {
					LUrow = _data[i];

					int maxel = Math.min(i, j);
					double s = 0.0;
					for (int k = 0; k < maxel; k++) {
						s += LUrow[k] * LUcol[k];
					}

					LUrow[j] = LUcol[i] -= s;
				}

				// Find pivot and swap if needed
				int p = j;
				for (int i = j + 1; i < rows; i++) {
					if (Double.compare(Math.abs(LUcol[i]), Math.abs(LUcol[p])) > 0) {
						p = i;
					}
				}
				if (p != j) {
					for (int k = 0; k < cols; k++) {
						double temp = _data[p][k];
						_data[p][k] = _data[j][k];
						_data[j][k] = temp;
					}
					int temp = _pivot[p];
					_pivot[p] = _pivot[j];
					_pivot[j] = temp;
				}
				// Wrapping up
				if (j < rows & Double.compare(_data[j][j], 0.0) != 0) {
					for (int i = j + 1; i < rows; i++) {
						_data[i][j] /= _data[j][j];
					}
				}
			}
		}

		public double det() {
			double det = _pivotsign;
			for (int i = 0; i < _cols; i++) {
				det *= _data[i][i];
			}
			return det;
		}

		public double[] solve(final double[] vector) {
			int rows = vector.length;
			double[] x = new double[rows];
			// Shuffle rows to match pivot vector
			for (int i = 0; i < rows; i++) {
				x[i] = vector[_pivot[i]];
			}

			// Solve Ly = b
			// Forward substitution
			for (int k = 0; k < rows; k++) {
				for (int i = k + 1; i < rows; i++) {
					x[i] -= x[k] * _data[i][k];
				}
			}

			// Solve Ux = y
			// Back substitution
			for (int k = rows - 1; k >= 0; k--) {
				x[k] /= _data[k][k];
				for (int i = 0; i < k; i++) {
					x[i] -= x[k] * _data[i][k];
				}
			}
			return x;
		}
	}

}
