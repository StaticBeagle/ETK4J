package com.wildbitsfoundry.etk4j.math.solvers.multivariate;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.function.Function;

import com.wildbitsfoundry.etk4j.constants.ETKConstants;
import com.wildbitsfoundry.etk4j.math.functions.MultivariateFunction;
import com.wildbitsfoundry.etk4j.util.NumArrays;

public class NewtonRaphson {
	
	public class SolverResults {
		double[] Solution;
		SolverStatus Status;
		int Iterations;
	}
	
	public enum SolverStatus {
		SUCCESS,
		DIVERGED,
		NOT_STARTED,
		MAX_VALUE_EXCEEDED,
		MIN_VALUE_EXCEEDED,
		ITERATION_LIMIT_EXCEEDED,
		
	}
	
	public enum ErrorEstimationScheme {
		SUM_ABS_ERROR,				// Norm 1
		MAX_ABS_ERROR,				// Infinite Norm
		SQRT_SUM_ABS_ERROR_SQUARED,	// Norm 2
	}
	
	private Solver _solver = null;

	public NewtonRaphson(MultivariateFunction[] functions, double[] initialGuess) {
		_solver = new RegularSolver(functions, initialGuess);
	}
	
	public NewtonRaphson(MultivariateFunction[] functions, MultivariateFunction[][] jacobian, double[] intialGuess) {
		_solver = new RegularSolver(functions, jacobian, intialGuess);
	}
	
	public NewtonRaphson(List<Function<Double[], Double>> functions, double[] initialGuess) {
		_solver = new FunctionalSolver(functions, initialGuess);
	}
	
	public NewtonRaphson(List<Function<Double[], Double>> functions, List<Function<Double[], Double[]>> jacobian, double[] intialGuess) {
		_solver = new FunctionalSolver(functions, jacobian, intialGuess);
	}

	public NewtonRaphson iterationLimit(int limit) {
		_solver.setMaxIter(limit);
		return this;
	}

	public NewtonRaphson absTolerance(double tol) {
		_solver.setAbsTol(tol);
		return this;
	}
	
	public NewtonRaphson relTolerance(double tol) {
		_solver.setRelTol(tol);
		return this;
	}

//	public NewtonRaphson maxAllowedValue(double val) {
//		_maxVal = val;
//		return this;
//	}

	public NewtonRaphson differentiationStepSize(double step) {
		_solver.setDiffStep(step);
		return this;
	}
	
	public NewtonRaphson setErrorEstimationScheme(ErrorEstimationScheme scheme) {
		_solver.setErrorEstimationScheme(scheme);
		return this;
	}
	
	public NewtonRaphson setConvergenceChecker(ConvergenceChecker checker) {
		_solver.setConvergenceChecker(checker);
		return this;
	}

	public double[] solve() {
		return _solver.solve();
	}
	
	private static void testRegular() {
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

		MultivariateFunction[][] jacobian = new MultivariateFunction[][] { { df1x1, df1x2, df1x3 },
				{ df2x1, df2x2, df2x3 }, { df3x1, df3x2, df3x3 } };

		MultivariateFunction[] functions = new MultivariateFunction[] { f1, f2, f3 };
		double[] initialguess = new double[] { 1d, 2d, 3d };

		System.out.printf("Solution with pre-computed Jacobian%n------------------------------------%n");
		double[] sol1 = new NewtonRaphson(functions, jacobian, initialguess)
				.absTolerance(1e-9)
				.relTolerance(1e-8)
				.iterationLimit(100)
				.solve();
		for (double val : sol1) {
			System.out.printf("%.6g%n", val);
		}

		System.out.printf("%nSolution with auto-computed Jacobian%n------------------------------------%n");
		double[] sol2 = new NewtonRaphson(functions, initialguess)
				.absTolerance(1e-9)
				.relTolerance(1e-8)
				.iterationLimit(100)
				.differentiationStepSize(0.001)
				.solve();
		for (double val : sol2) {
			System.out.printf("%.6g%n", val);
		}
	}
	
	public static void testFunctional() {
		double[] initialguess = { 1d, 2d, 3d };

		List<Function<Double[], Double>> functions = new ArrayList<>();
		functions.add(x -> 3 * x[0] - Math.cos(x[1] * x[2]) - 1.5);
		functions.add(x -> 4 * x[0] * x[0] - 625.0 * x[1] * x[1] + 2.0 * x[2] - 1.0);
		functions.add(x -> 20 * x[2] + Math.exp(-x[0] * x[1]) + 9.0);

		List<Function<Double[], Double[]>> jacobian = new ArrayList<>();
		jacobian.add(x -> new Double[] { 3.0, x[2] * Math.sin(x[1] * x[2]), x[1] * Math.sin(x[1] * x[2]) });
		jacobian.add(x -> new Double[] { 8 * x[0], -1250.0 * x[1], 2.0 });
		jacobian.add(x -> new Double[] { -x[1] * Math.exp(x[0] * x[1]), -x[0] * Math.exp(x[0] * x[1]), 20.0 });

		System.out.printf("Solution with pre-computed Jacobian%n------------------------------------%n");
		double[] sol1 = new NewtonRaphson(functions, jacobian, initialguess)
				.absTolerance(1e-9)
				.relTolerance(1e-6)
				.iterationLimit(100)
				.solve();
		for (double val : sol1) {
			System.out.printf("%.6g%n", val);
		}

		System.out.printf("%nSolution with auto-computed Jacobian%n------------------------------------%n");
		double[] sol2 = new NewtonRaphson(functions, initialguess)
				.absTolerance(1e-9)
				.relTolerance(1e-6)
				.iterationLimit(100)
				.differentiationStepSize(0.001)
				.solve();
		for (double val : sol2) {
			System.out.printf("%.6g%n", val);
		}
	}
	
	public static void main(String[] args) {
		testRegular();
		testFunctional();
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
		
		public void solve(final double[] vector, double[] solution) {
			int rows = vector.length;
			// Shuffle rows to match pivot vector
			for (int i = 0; i < rows; i++) {
				solution[i] = vector[_pivot[i]];
			}

			// Solve Ly = b
			// Forward substitution
			for (int k = 0; k < rows; k++) {
				for (int i = k + 1; i < rows; i++) {
					solution[i] -= solution[k] * _data[i][k];
				}
			}

			// Solve Ux = y
			// Back substitution
			for (int k = rows - 1; k >= 0; k--) {
				solution[k] /= _data[k][k];
				for (int i = 0; i < k; i++) {
					solution[i] -= solution[k] * _data[i][k];
				}
			}
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
	
	private interface Solver {
		public double[] solve();
		public void setMaxIter(int max);
		public void setAbsTol(double tol);
		public void setRelTol(double tol);
		public void setDiffStep(double step);
		public void setConvergenceChecker(ConvergenceChecker checker);
		public void setErrorEstimationScheme(ErrorEstimationScheme scheme);
	}
	
	private interface ConvergenceChecker {
		public boolean checkForConvergence(double[] x, double[] y, double absTol, double relTol);
		
		default public boolean checkForConvergence(Double[] x, Double[] y, double absTol, double relTol) {
			return checkForConvergence(NumArrays.unbox(x), NumArrays.unbox(y), absTol, relTol);
		}
	}
	
	private static abstract class NewtonRaphsonSolver <T, U> implements Solver {
		protected int _maxIter = 100;
		protected double _absTol = 1e-9;
		protected double _relTol = 1e-6;
		protected double _maxVal = Double.POSITIVE_INFINITY;
		protected double _step = 0.001;
		protected double[] _x0;
		
		protected T _functions;
		protected U _jacobian;
		private ConvergenceChecker _convChecker;
		
		protected NewtonRaphsonSolver(T functions, double[] x0) {
			this(functions, null, x0);
		}
		
		protected NewtonRaphsonSolver(T functions, U jacobian, double[] x0) {
			_functions = functions;
			_jacobian = jacobian;
			_x0 = x0;
			this.setErrorEstimationScheme(ErrorEstimationScheme.SQRT_SUM_ABS_ERROR_SQUARED);
		}
		
		public abstract double[] solve(Optional<U> jacobian);
		
		@Override
		public void setMaxIter(int max) {
			_maxIter = max;
		}
		
		@Override
		public void setAbsTol(double tol) {
			_absTol = tol;
		}
		
		@Override
		public void setRelTol(double tol) {
			_relTol = tol;
		}
		
		@Override
		public void setDiffStep(double step) {
			_step = step;
		}
		
		@Override
		public final double[] solve() {
			return this.solve(Optional.ofNullable(_jacobian));
		}
		
		@Override
		public void setErrorEstimationScheme(ErrorEstimationScheme scheme) {
			switch(scheme) {
			case MAX_ABS_ERROR:
				_convChecker = new ConvergenceChecker() {
					
					@Override
					public boolean checkForConvergence(double[] x, double[] y, double absTol, double relTol) {
						double delta = NumArrays.normInf(NumArrays.subtract(x, y));
						double tmp = Math.min(NumArrays.normInf(x), NumArrays.normInf(y));
						return delta < absTol + relTol * tmp;
					}
				};
				break;
			case SQRT_SUM_ABS_ERROR_SQUARED:
				_convChecker = new ConvergenceChecker() {
					
					@Override
					public boolean checkForConvergence(double[] x, double[] y, double absTol, double relTol) {
						double delta = NumArrays.distance(x, y);
						double tmp = Math.min(NumArrays.norm(x), NumArrays.norm(y));
						return delta < absTol + relTol * tmp;
					}
				};
				break;
			case SUM_ABS_ERROR:
				_convChecker = new ConvergenceChecker() {
					
					@Override
					public boolean checkForConvergence(double[] x, double[] y, double absTol, double relTol) {
						double delta = NumArrays.norm1(NumArrays.subtract(x, y));
						double tmp = Math.min(NumArrays.norm1(x), NumArrays.norm1(y));
						return delta < absTol + relTol * tmp;
					}
				};
				break;
			default:
				break;
			}
		}
		
		@Override
		public void setConvergenceChecker(ConvergenceChecker checker) {
			_convChecker = checker;
		}
		
		protected boolean checkForConvergence(double[] xfinal, double[] xcurrent) {
			return _convChecker.checkForConvergence(xfinal, xcurrent, _absTol, _relTol);
		}
		
		protected boolean checkForConvergence(Double[] xfinal, Double[] xcurrent) {
			return _convChecker.checkForConvergence(xfinal, xcurrent, _absTol, _relTol);
		}
	}
	
	private static class RegularSolver extends NewtonRaphsonSolver<MultivariateFunction[], MultivariateFunction[][]> {

		public RegularSolver(MultivariateFunction[] functions, double[] x0) {
			super(functions, x0);
		}
		
		protected RegularSolver(MultivariateFunction[] functions, MultivariateFunction[][] jacobian, double[] x0) {
			super(functions, jacobian, x0);
		}
		
		private static void evaluateJacobian(double[] x, MultivariateFunction[][] jacobian, double[][] out) {
			int rows = out.length;
			int cols = out[0].length;
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					out[i][j] = jacobian[i][j].evaluateAt(x);
				}
			}
		}

		private static void evaluateFunctionsAt(double[] x, MultivariateFunction[] functions, double[] out) {
			int length = x.length;
			for (int i = 0; i < length; i++) {
				out[i] = functions[i].evaluateAt(x);
			}
		}
		
		private static void computeJacobian(double[] x, MultivariateFunction[] f, double[] dfx, double[] fx, double step, double[][] out) {
			double h = 1.0 / step;
			int rows = f.length;
			int cols = x.length;
			for (int i = 0; i < cols; i++) {
				x[i] = x[i] + step;
				evaluateFunctionsAt(x, f, dfx);
				x[i] = x[i] - step;
				evaluateFunctionsAt(x, f, fx);
				for (int j = 0; j < rows; j++) {
					out[j][i] = (dfx[j] - fx[j]) * h;
				}
			}
		}

		@Override
		public double[] solve(Optional<MultivariateFunction[][]> jacobian) {
			int maxiter = _maxIter;
			double maxval = _maxVal;
			double step = _step;
			double[] x0 = _x0;
			MultivariateFunction[] functions = _functions;

			int rows = functions.length;
			int cols = x0.length;
			double[] xcurrent = new double[cols];
			double[] xfinal = new double[cols];
			double[] residuals = new double[cols];
			double[] functionValues = new double[cols];
			double[][] jacobianMatrixValues = new double[rows][cols];

			System.arraycopy(x0, 0, xcurrent, 0, rows);
			while (maxiter-- > 0) {
				if(jacobian.isPresent()) {
					evaluateJacobian(xcurrent, jacobian.get(), jacobianMatrixValues);
					evaluateFunctionsAt(xcurrent, functions, functionValues);
				} else {
					computeJacobian(xcurrent, functions, residuals, functionValues, step, jacobianMatrixValues);				
				}
				
				LU matrix = new LU(jacobianMatrixValues);
				if (Double.compare(Math.abs(matrix.det()), ETKConstants.DOUBLE_EPS) < 0) {
					// bail out the matrix is singular
					return null; // <-- we'll think about this later
				}

				matrix.solve(functionValues, residuals);
				for (int i = 0; i < cols; i++) {
					xfinal[i] = xcurrent[i] - residuals[i];
				}
				evaluateFunctionsAt(xfinal, functions, functionValues);

				if(this.checkForConvergence(xfinal, xcurrent)) {
					return xfinal;
				}

				boolean diverged = false;
				for (int i = 0; i < cols; i++) {
					diverged |= Double.compare(Math.abs(xcurrent[i]), maxval) > 0;
				}
				if (diverged) {
					return null;
				}
				System.arraycopy(xfinal, 0, xcurrent, 0, rows);
			}
			return null;
		}		
	}

	private static class FunctionalSolver extends NewtonRaphsonSolver<List<Function<Double[], Double>>, List<Function<Double[], Double[]>>> {

		protected FunctionalSolver(List<Function<Double[], Double>> functions, double[] x0) {
			super(functions, x0);
		}
		
		protected FunctionalSolver(List<Function<Double[], Double>> functions, List<Function<Double[], Double[]>> jacobian, double[] x0) {
			super(functions, jacobian, x0);
		}
		
		private static void computeJacobian(Double[] x, List<Function<Double[], Double>> f, double[] dfx, double[] fx, double step, double[][] out) {
			double h = 1.0 / step;
			int rows = f.size();
			int cols = x.length;
			for (int i = 0; i < cols; i++) {
				x[i] = x[i] + step;
				evaluateFunctionsAt(x, f, dfx);
				x[i] = x[i] - step;
				evaluateFunctionsAt(x, f, fx);
				for (int j = 0; j < rows; j++) {
					out[j][i] = (dfx[j] - fx[j]) * h;
				}
			}
		}
		
		private static void evaluateJacobian(Double[] xcurrent, List<Function<Double[], Double[]>> jacobian, double[][] out) {
			int rows = out.length;

			for (int i = 0; i < rows; ++i) {
				Double[] tmp = jacobian.get(i).apply(xcurrent);
				for(int j = 0; j < tmp.length; ++j){
					out[i][j] = tmp[j].doubleValue();
				}
			}
		}

		private static void evaluateFunctionsAt(Double[] x, List<Function<Double[], Double>> functions, double[] out) {
			int length = x.length;
			for (int i = 0; i < length; i++) {
				out[i] = functions.get(i).apply(x).doubleValue();
			}
		}

		@Override
		public double[] solve(Optional<List<Function<Double[], Double[]>>> jacobian) {
			int maxiter = _maxIter;
			double tol = _absTol;
			double maxval = _maxVal;
			double step = _step;
			Double[] x0 = NumArrays.box(_x0);
			List<Function<Double[], Double>> functions = _functions;

			int rows = functions.size();
			int cols = x0.length;
			Double[] xcurrent = new Double[cols];
			Double[] xfinal = new Double[cols];
			double[] residuals = new double[cols];
			double[] functionValues = new double[cols];
			double[][] jacobianMatrixValues = new double[rows][cols];

			System.arraycopy(x0, 0, xcurrent, 0, rows);
			while (maxiter-- > 0) {

				if (jacobian.isPresent()) {
					evaluateJacobian(xcurrent, jacobian.get(), jacobianMatrixValues);
					evaluateFunctionsAt(xcurrent, functions, functionValues);
				} else {
					 computeJacobian(xcurrent, functions, residuals, functionValues, step, jacobianMatrixValues);
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
				evaluateFunctionsAt(xfinal, functions ,functionValues);

				if (this.checkForConvergence(xfinal, xcurrent)) {
					return NumArrays.unbox(xfinal);
				}

				boolean diverged = false;
				for (int i = 0; i < cols; i++) {
					diverged |= Double.compare(Math.abs(xcurrent[i]), maxval) > 0;
				}
				if (diverged) {
					return null;
				}

				System.arraycopy(xfinal, 0, xcurrent, 0, rows);
			}
			return null;
		}
	}
}