package com.wildbitsfoundry.etk4j.math.solvers.univariate;

public interface UnivariateSolver {
	public SolverResults solve();
	public void setMaxIter(int max);
	public void setAbsTol(double tol);
	public void setRelTol(double tol);
	public void setMaxAbsVal(double max);
	public void setDiffStep(double step);
}
