package com.wildbitsfoundry.etk4j.math.solvers.multivariate;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;

public class NewtonRaphsonFunctional {

	private int _maxIter = 100;
	private double _absTol = 1e-9;
	private double _maxVal = Double.POSITIVE_INFINITY;
	private boolean _autoJacobian = true;
	private double _step = 0.001;
	private Double[] _x0;

	private List<Function<Double[], Double>> _functions;
	private List<Function<Double[], Double[]>> _jacobian;

	public NewtonRaphsonFunctional(List<Function<Double[], Double>> functions) {
		_functions = functions;
	}

	public NewtonRaphsonFunctional(List<Function<Double[], Double>> functions,
			List<Function<Double[], Double[]>> jacobian) {
		_functions = functions;
		_jacobian = jacobian;
		_autoJacobian = false;
	}

	public NewtonRaphsonFunctional iterationLimit(int limit) {
		_maxIter = limit;
		return this;
	}

	public NewtonRaphsonFunctional absTolerance(double tol) {
		_absTol = tol;
		return this;
	}

	public NewtonRaphsonFunctional maxAllowedValue(double val) {
		_maxVal = val;
		return this;
	}

	public NewtonRaphsonFunctional differentiationStepSize(double step) {
		_step = step;
		return this;
	}

	public NewtonRaphsonFunctional initialGuess(Double[] x0) {
		_x0 = Arrays.copyOf(x0, x0.length);
		return this;
	}

	private Double[] solvePreComputedJacobian() {
		int maxiter = _maxIter;
		double tol = _absTol;
		double maxval = _maxVal;
		Double[] x0 = _x0;
		List<Function<Double[], Double>> functions = _functions;
		List<Function<Double[], Double[]>> jacobian = _jacobian;

		int rows = jacobian.size();
		int cols = jacobian.size();
		Double[] xcurrent = x0;
		Double[] xfinal = new Double[cols];
		double[] functionValues = new double[cols];
		Double[][] jacobianMatrixValues = new Double[rows][cols];
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

	public Double[] solve() {
		return _autoJacobian ? this.solveAutoComputeJacobian() : this.solvePreComputedJacobian();
	}

	private Double[] solveAutoComputeJacobian() {
		int maxiter = _maxIter;
		double tol = _absTol;
		double maxval = _maxVal;
		double step = _step;
		double h = 1.0 / step;
		Double[] x0 = _x0;
		List<Function<Double[], Double>> functions = _functions;

		int rows = functions.size();
		int cols = x0.length;
		Double[] xcurrent = new Double[cols];
		Double[] xfinal = new Double[cols];
		double[] residuals = new double[cols];
		double[] functionValues = new double[cols];
		double[] functionPrimeValues = new double[cols];
		Double[][] jacobianMatrixValues = new Double[rows][cols];

		System.arraycopy(x0, 0, xcurrent, 0, rows);
		while (maxiter > 0) {
			// Compute Jacobian
			for (int i = 0; i < cols; i++) {
				xcurrent[i] = xcurrent[i] + step;
				eval(xcurrent, functionPrimeValues, functions);
				xcurrent[i] = xcurrent[i] - step;
				eval(xcurrent, functionValues, functions);
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
				xfinal[i] = xcurrent[i] - residuals[i];
			}
			eval(xfinal, functionValues, functions);

			double relErrorNorm = 0.0;
			double errorNorm = 0.0;
			for (int i = 0; i < cols; i++) {
				// calculate the norm2 of the relative error
				relErrorNorm += (xfinal[i] - xcurrent[i]) * (xfinal[i] - xcurrent[i]);
				errorNorm += xfinal[i] * xfinal[i];
			}
			double error = Math.sqrt(relErrorNorm) / Math.sqrt(errorNorm);
			if (Double.compare(error, tol) <= 0) {
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
			System.arraycopy(xfinal, 0, xcurrent, 0, rows);
		}
		return null;
	}

	private static void evaluateJacobian(Double[] xcurrent, Double[][] jacobianMatrixValues,
			List<Function<Double[], Double[]>> jacobian) {
		int rows = jacobianMatrixValues.length;

		for (int i = 0; i < rows; i++) {
			jacobianMatrixValues[i] = jacobian.get(i).apply(xcurrent);

		}
	}

	private static void eval(Double[] x, double[] out, List<Function<Double[], Double>> functions) {
		int length = x.length;
		for (int i = 0; i < length; i++) {
			out[i] = functions.get(i).apply(x);
		}
	}

	public static void main(String[] args) {

		Double[] initialguess = { 1d, 2d, 3d };

		List<Function<Double[], Double>> functions = new ArrayList<>();
		functions.add(x -> 3 * x[0] - Math.cos(x[1] * x[2]) - 1.5);
		functions.add(x -> 4 * x[0] * x[0] - 625.0 * x[1] * x[1] + 2.0 * x[2] - 1.0);
		functions.add(x -> 20 * x[2] + Math.exp(-x[0] * x[1]) + 9.0);

		List<Function<Double[], Double[]>> jacobian = new ArrayList<>();
		jacobian.add(x -> new Double[] { 3.0, x[2] * Math.sin(x[1] * x[2]), x[1] * Math.sin(x[1] * x[2]) });
		jacobian.add(x -> new Double[] { 8 * x[0], -1250.0 * x[1], 2.0 });
		jacobian.add(x -> new Double[] { -x[1] * Math.exp(x[0] * x[1]), -x[0] * Math.exp(x[0] * x[1]), 20.0 });

		System.out.printf("Solution with pre-computed Jacobian%n------------------------------------%n");
		Double[] sol1 = new NewtonRaphsonFunctional(functions, jacobian).absTolerance(1e-9).iterationLimit(100)
				.initialGuess(initialguess).solve();
		for (double val : sol1) {
			System.out.printf("%.6g%n", val);
		}

		System.out.printf("%nSolution with auto-computed Jacobian%n------------------------------------%n");
		Double[] sol2 = new NewtonRaphsonFunctional(functions).absTolerance(1e-9).iterationLimit(100)
				.initialGuess(initialguess).differentiationStepSize(0.001).solve();
		for (double val : sol2) {
			System.out.printf("%.6g%n", val);
		}
	}

	private static class LU {
		protected Double[][] _data;
		protected int _cols;

		protected int _pivotsign;
		protected int[] _pivot;

		public LU(Double[][] matrix) {
			int rows = matrix.length;
			int cols = _cols = matrix[0].length;
			_data = matrix;

			_pivot = new int[rows];
			for (int i = 0; i < rows; i++) {
				_pivot[i] = i;
			}
			_pivotsign = 1;

			Double[] LUrow;
			Double[] LUcol = new Double[rows];

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
