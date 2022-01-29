package com.wildbitsfoundry.etk4j.math.optimize.solvers;

public class SolverResults<T> {
	private T value;
	private int numberOfIterations;
	private String solverStatus;
	private double error;
	private boolean converged;

	public T getValue() {
		return value;
	}

	void setValue(T value) {
		this.value = value;
	}

	public int getNumberOfIterations() {
		return numberOfIterations;
	}

	void setNumberOfIterations(int numberOfIterations) {
		this.numberOfIterations = numberOfIterations;
	}

	public String getSolverStatus() {
		return solverStatus;
	}

	void setSolverStatus(String solverStatus) {
		this.solverStatus = solverStatus;
	}

	public double getError() {
		return error;
	}

	void setError(double error) {
		this.error = error;
	}

	public boolean hasConverged() {
		return converged;
	}

	void setHasConverged(boolean converged) {
		this.converged = converged;
	}

	@Override
	public String toString() {
		return "SolverResults{" +
				"value=" + value +
				", numberOfIterations=" + numberOfIterations +
				", solverStatus='" + solverStatus + '\'' +
				", error=" + error +
				", converged=" + converged +
				'}';
	}
}

