package com.wildbitsfoundry.etk4j.math.solvers.multivariate;

import com.wildbitsfoundry.etk4j.math.solvers.multivariate.ErrorEstimationScheme;

public interface MultivariateSolver {
	public SolverResults solve();
	public MultivariateSolver setAbsTolerance(double tol);
	public MultivariateSolver setRelTolerance(double tol);
	public MultivariateSolver setIterationLimit(int limit);
	public MultivariateSolver setMaxAbsAllowedValue(double max);
	public MultivariateSolver setDifferentiationStepSize(double step);
	public MultivariateSolver setErrorEstimationScheme(ErrorEstimationScheme scheme);
}


