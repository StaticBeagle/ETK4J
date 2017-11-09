package com.wildbitsfoundry.etk4j.math.solvers.multivariate;

import com.wildbitsfoundry.etk4j.math.solvers.multivariate.ErrorEstimationScheme;

public interface MultivariateSolver {
	public SolverResults solve();
	public void setMaxIter(int max);
	public void setAbsTol(double tol);
	public void setRelTol(double tol);
	public void setMaxAbsVal(double max);
	public void setDiffStep(double step);
	public void setConvergenceChecker(ConvergenceChecker checker);
	public void setErrorEstimationScheme(ErrorEstimationScheme scheme);
}
