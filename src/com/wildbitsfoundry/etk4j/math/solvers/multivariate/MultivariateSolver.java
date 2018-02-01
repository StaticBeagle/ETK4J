package com.wildbitsfoundry.etk4j.math.solvers.multivariate;

import com.wildbitsfoundry.etk4j.math.solvers.multivariate.ErrorEstimationScheme;

public interface MultivariateSolver {
	public SolverResults solve();
	public MultivariateSolver setMaxIter(int max);
	public MultivariateSolver setAbsTol(double tol);
	public MultivariateSolver setRelTol(double tol);
	public MultivariateSolver setMaxAbsVal(double max);
	public MultivariateSolver setDiffStep(double step);
	public MultivariateSolver setConvergenceChecker(ConvergenceChecker checker);
	public MultivariateSolver setErrorEstimationScheme(ErrorEstimationScheme scheme);
}
