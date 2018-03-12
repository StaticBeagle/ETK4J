package com.wildbitsfoundry.etk4j.math.solvers.multivariate;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.functions.MultivariateFunction;
import com.wildbitsfoundry.etk4j.math.solvers.multivariate.SolverResults.SolverStatus;
import com.wildbitsfoundry.etk4j.math.solvers.multivariate.SolverUtils.LU;
import com.wildbitsfoundry.etk4j.util.NumArrays;

public class NewtonRaphson implements MultivariateSolver {
	
	protected int _maxIter = 100;
	protected double _absTol = 1e-9;
	protected double _relTol = 1e-6;
	protected double _maxVal = Double.NaN;
	protected double _step = 1e-7;
	protected double[] _x0;
	
	MultivariateFunction[] _functions;
	MultivariateFunction[][] _jacobian;
	
	private ErrorEstimationScheme _errorScheme;


	public NewtonRaphson(MultivariateFunction[] functions, double[] initialGuess) {
		this(functions, null, initialGuess);
	}

	public NewtonRaphson(MultivariateFunction[] functions, MultivariateFunction[][] jacobian, double[] initialGuess) {
		_functions = functions;
		_jacobian = jacobian;
		_x0 = initialGuess;
		_errorScheme = ErrorEstimationScheme.SQRT_SUM_ABS_ERROR_SQUARED;
	}

	@Override
	public SolverResults solve() {
		return this.solve(_jacobian);
	}

	@Override
	public NewtonRaphson setIterationLimit(int limit) {
		_maxIter = limit;
		return this;
	}

	@Override
	public NewtonRaphson setAbsTolerance(double tol) {
		_absTol = tol;
		return this;
	}

	@Override
	public NewtonRaphson setRelTolerance(double tol) {
		_relTol = tol;
		return this;
	}

	@Override
	public NewtonRaphson setMaxAbsAllowedValue(double max) {
		_maxVal = max;
		return this;
	}

	@Override
	public NewtonRaphson setDifferentiationStepSize(double step) {
		_step = step;
		return this;
	}

	@Override
	public NewtonRaphson setErrorEstimationScheme(ErrorEstimationScheme scheme) {
		_errorScheme = scheme;
		return this;
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

	protected SolverResults solve(MultivariateFunction[][] jacobian) {
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
			if (jacobian != null) {
				evaluateJacobian(xcurrent, jacobian, jacobianMatrixValues);
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

	
	protected double estimateError(double[] x, double[] y) {
		return _errorScheme.calculateRelativeError(x, y);
	}

	protected boolean checkForConvergence(double[] x, double[] y) {
		double delta = this.estimateError(x, y);
		double tmp = _errorScheme.calculateMaxNormOfError(x, y);
		return delta < _absTol + _relTol * tmp;
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
}