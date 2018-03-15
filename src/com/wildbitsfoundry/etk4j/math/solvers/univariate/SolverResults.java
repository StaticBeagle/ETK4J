package com.wildbitsfoundry.etk4j.math.solvers.univariate;


public class SolverResults {
	public enum SolverStatus {
		SUCCESS,
		DIVERGED,
		NOT_STARTED,
		JACOBIAN_IS_SINGULAR,
		MAX_ABS_VALUE_EXCEEDED,
		ITERATION_LIMIT_EXCEEDED,
	}
	
	public double Solution;
	public SolverStatus Status;
	public int Iterations;
	boolean Converged;
	double EstimatedError;
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(String.format("Status: %s%n", Status))
			.append(String.format("Iterations: %d%n", Iterations))
			.append(String.format("Converged: %s%n", Converged))
			.append(String.format("Estimated Error: %s%n", EstimatedError))
			.append(String.format("Solution: %.6g", Solution));

		return sb.toString();
	}
	
}

