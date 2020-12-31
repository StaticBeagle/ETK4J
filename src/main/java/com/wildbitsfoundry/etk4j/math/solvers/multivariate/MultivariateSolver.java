package com.wildbitsfoundry.etk4j.math.solvers.multivariate;

import com.wildbitsfoundry.etk4j.math.solvers.multivariate.ErrorEstimationScheme;

public interface MultivariateSolver {
	SolverResults solve();
	MultivariateSolver setAbsTolerance(double tol);
	MultivariateSolver setRelTolerance(double tol);
	MultivariateSolver setIterationLimit(int limit);
	MultivariateSolver setMaxAbsAllowedValue(double max);
	MultivariateSolver setDifferentiationStepSize(double step);
	MultivariateSolver setErrorEstimationScheme(ErrorEstimationScheme scheme);
}


