package com.wildbitsfoundry.etk4j.math.solvers.multivariate;

import java.util.List;
import java.util.Optional;
import java.util.function.Function;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.functions.MultivariateFunction;
import com.wildbitsfoundry.etk4j.math.solvers.multivariate.SolverResults.SolverStatus;
import com.wildbitsfoundry.etk4j.util.NumArrays;

import static com.wildbitsfoundry.etk4j.math.solvers.multivariate.SolverUtils.LU;

public class NewtonRaphson {

	private MultivariateSolver _solver = null;

	public NewtonRaphson(MultivariateFunction[] functions, double[] initialGuess) {
		_solver = new RegularSolver(functions, initialGuess);
	}

	public NewtonRaphson(MultivariateFunction[] functions, MultivariateFunction[][] jacobian, double[] intialGuess) {
		_solver = new RegularSolver(functions, jacobian, intialGuess);
	}

	public NewtonRaphson(List<Function<Double[], Double>> functions, double[] initialGuess) {
		_solver = new FunctionalSolver(functions, initialGuess);
	}

	public NewtonRaphson(List<Function<Double[], Double>> functions, List<Function<Double[], Double[]>> jacobian,
			double[] intialGuess) {
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

	public NewtonRaphson maxAbsAllowedValue(double max) {
		_solver.setMaxAbsVal(max);
		return this;
	}

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

	public SolverResults solve() {
		return _solver.solve();
	}

	private static abstract class NewtonRaphsonSolver<T, U> implements MultivariateSolver {
		protected int _maxIter = 100;
		protected double _absTol = 1e-9;
		protected double _relTol = 1e-6;
		protected double _maxVal = Double.NaN;
		protected double _step = 1e-7;
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

		public abstract SolverResults solve(Optional<U> jacobian);

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
		public final SolverResults solve() {
			return this.solve(Optional.ofNullable(_jacobian));
		}

		@Override
		public void setMaxAbsVal(double max) {
			_maxVal = max;
		}

		@Override
		public void setErrorEstimationScheme(ErrorEstimationScheme scheme) {
			switch (scheme) {
			case MAX_ABS_ERROR:
				_convChecker = new ConvergenceChecker() {
					@Override
					public double estimateError(double[] x, double[] y) {
						return NumArrays.normInf(NumArrays.subtract(x, y));
					}

					@Override
					public boolean checkForConvergence(double[] x, double[] y, double absTol, double relTol) {
						double delta = this.estimateError(x, y);
						double tmp = Math.min(NumArrays.normInf(x), NumArrays.normInf(y));
						return delta < absTol + relTol * tmp;
					}
				};
				break;
			case SQRT_SUM_ABS_ERROR_SQUARED:
				_convChecker = new ConvergenceChecker() {
					@Override
					public double estimateError(double[] x, double[] y) {
						return NumArrays.distance(x, y);
					}

					@Override
					public boolean checkForConvergence(double[] x, double[] y, double absTol, double relTol) {
						double delta = this.estimateError(x, y);
						double tmp = Math.min(NumArrays.norm2(x), NumArrays.norm2(y));
						return delta < absTol + relTol * tmp;
					}
				};
				break;
			case SUM_ABS_ERROR:
				_convChecker = new ConvergenceChecker() {
					@Override
					public double estimateError(double[] x, double[] y) {
						return NumArrays.norm1(NumArrays.subtract(x, y));
					}

					@Override
					public boolean checkForConvergence(double[] x, double[] y, double absTol, double relTol) {
						double delta = this.estimateError(x, y);
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
		
		protected double estimateError(double[] xfinal, double[] xcurrent) {
			return _convChecker.estimateError(xfinal, xcurrent);
		}

		protected boolean checkForConvergence(double[] xfinal, double[] xcurrent) {
			return _convChecker.checkForConvergence(xfinal, xcurrent, _absTol, _relTol);
		}

		protected boolean checkForConvergence(Double[] xfinal, Double[] xcurrent) {
			return _convChecker.checkForConvergence(xfinal, xcurrent, _absTol, _relTol);
		}

		protected static SolverResults buildResults(double[] xfinal, SolverStatus status, int iterCount, double error) {
			SolverResults sr = new SolverResults();
			sr.Iterations = iterCount;
			sr.Solution = xfinal;
			sr.Status = status;
			sr.Converged = status.equals(SolverStatus.SUCCESS) ? true : false;
			sr.EstimatedError = error;
			return sr;
		}

		protected static SolverResults buildResults(Double[] xfinal, SolverStatus status, int iterCount, double error) {
			return buildResults(NumArrays.unbox(xfinal), status, iterCount, error);
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

		private static void computeJacobian(double[] x, MultivariateFunction[] f, double[] dfx, double[] fx,
				double step, double[][] out) {
			double den = 0.5 / step;
			int rows = f.length;
			int cols = x.length;
			for (int i = 0; i < cols; i++) {
				x[i] = x[i] + step;
				evaluateFunctionsAt(x, f, dfx);
				x[i] = x[i] - 2.0 * step;
				evaluateFunctionsAt(x, f, fx);
				for (int j = 0; j < rows; j++) {
					out[j][i] = (dfx[j] - fx[j]) * den;
				}
			}
		}

		@Override
		public SolverResults solve(Optional<MultivariateFunction[][]> jacobian) {
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
				if (jacobian.isPresent()) {
					evaluateJacobian(xcurrent, jacobian.get(), jacobianMatrixValues);
					evaluateFunctionsAt(xcurrent, functions, functionValues);
				} else {
					computeJacobian(xcurrent, functions, residuals, functionValues, step, jacobianMatrixValues);
				}

				LU matrix = new LU(jacobianMatrixValues);
				if (Double.compare(Math.abs(matrix.det()), ConstantsETK.DOUBLE_EPS) <= 0) {
					double error = this.estimateError(xfinal, xcurrent);
					return buildResults(xfinal, SolverStatus.JACOBIAN_IS_SINGULAR, _maxIter - maxiter, error);
				}

				matrix.solve(functionValues, residuals);
				for (int i = 0; i < cols; i++) {
					xfinal[i] = xcurrent[i] - residuals[i];
				}
				evaluateFunctionsAt(xfinal, functions, functionValues);

				if (this.checkForConvergence(xfinal, xcurrent)) {
					double error = this.estimateError(xfinal, xcurrent);
					return buildResults(xfinal, SolverStatus.SUCCESS, _maxIter - maxiter, error);
				}

				if (!Double.isNaN(maxval) && Double.compare(NumArrays.normInf(xfinal), maxval) >= 0) {
					double error = this.estimateError(xfinal, xcurrent);
					return buildResults(xfinal, SolverStatus.MAX_ABS_VALUE_EXCEEDED, _maxIter - maxiter, error);
				}
				System.arraycopy(xfinal, 0, xcurrent, 0, rows);
			}
			double error = this.estimateError(xfinal, xcurrent);
			return buildResults(xfinal, SolverStatus.ITERATION_LIMIT_EXCEEDED, _maxIter - maxiter, error);
		}
	}

	private static class FunctionalSolver
			extends NewtonRaphsonSolver<List<Function<Double[], Double>>, List<Function<Double[], Double[]>>> {

		protected FunctionalSolver(List<Function<Double[], Double>> functions, double[] x0) {
			super(functions, x0);
		}

		protected FunctionalSolver(List<Function<Double[], Double>> functions,
				List<Function<Double[], Double[]>> jacobian, double[] x0) {
			super(functions, jacobian, x0);
		}

		private static void computeJacobian(Double[] x, List<Function<Double[], Double>> f, double[] dfx, double[] fx,
				double step, double[][] out) {
			double h = 0.5 / step;
			int rows = f.size();
			int cols = x.length;
			for (int i = 0; i < cols; i++) {
				x[i] = x[i] + step;
				evaluateFunctionsAt(x, f, dfx);
				x[i] = x[i] - 2.0 * step;
				evaluateFunctionsAt(x, f, fx);
				for (int j = 0; j < rows; j++) {
					out[j][i] = (dfx[j] - fx[j]) * h;
				}
			}
		}

		private static void evaluateJacobian(Double[] xcurrent, List<Function<Double[], Double[]>> jacobian,
				double[][] out) {
			int rows = out.length;

			for (int i = 0; i < rows; ++i) {
				Double[] tmp = jacobian.get(i).apply(xcurrent);
				for (int j = 0; j < tmp.length; ++j) {
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
		public SolverResults solve(Optional<List<Function<Double[], Double[]>>> jacobian) {
			int maxiter = _maxIter;
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
				if (Double.compare(Math.abs(matrix.det()), ConstantsETK.DOUBLE_EPS) <= 0) {
					double error = this.estimateError(NumArrays.unbox(xfinal), NumArrays.unbox(xcurrent));
					return buildResults(xfinal, SolverStatus.JACOBIAN_IS_SINGULAR, _maxIter - maxiter, error);
				}

				residuals = matrix.solve(functionValues);
				for (int i = 0; i < cols; i++) {
					xfinal[i] = xcurrent[i] - residuals[i];
				}
				evaluateFunctionsAt(xfinal, functions, functionValues);

				if (this.checkForConvergence(xfinal, xcurrent)) {
					double error = this.estimateError(NumArrays.unbox(xfinal), NumArrays.unbox(xcurrent));
					return buildResults(xfinal, SolverStatus.SUCCESS, _maxIter - maxiter, error);
				}

				if (!Double.isNaN(maxval) && Double.compare(NumArrays.normInf(xfinal), maxval) >= 0) {
					double error = this.estimateError(NumArrays.unbox(xfinal), NumArrays.unbox(xcurrent));
					return buildResults(xfinal, SolverStatus.MAX_ABS_VALUE_EXCEEDED, _maxIter - maxiter, error);
				}

				System.arraycopy(xfinal, 0, xcurrent, 0, rows);
			}
			double error = this.estimateError(NumArrays.unbox(xfinal), NumArrays.unbox(xcurrent));
			return buildResults(xfinal, SolverStatus.ITERATION_LIMIT_EXCEEDED, _maxIter - maxiter, error);
		}
	}
}