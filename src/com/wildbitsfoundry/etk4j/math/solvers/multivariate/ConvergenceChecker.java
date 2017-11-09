package com.wildbitsfoundry.etk4j.math.solvers.multivariate;

import com.wildbitsfoundry.etk4j.util.NumArrays;

public interface ConvergenceChecker {
	public double estimateError(double[] x, double[] y);
	public boolean checkForConvergence(double[] x, double[] y, double absTol, double relTol);
	
	default public boolean checkForConvergence(Double[] x, Double[] y, double absTol, double relTol) {
		return checkForConvergence(NumArrays.unbox(x), NumArrays.unbox(y), absTol, relTol);
	}
}