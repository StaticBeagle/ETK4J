package com.wildbitsfoundry.etk4j.math.optimize;


public class SolverResults<T> {
	private T value;
	private int numberOfIterations;
	private SolverStatus solverStatus;
	private double error;

	public enum SolverStatus {
		CONVERGED("Converged"),
		MAX_VALUE_EXCEEDED("Maximum allowed value exceeded"),
		MIN_VALUE_EXCEEDED("Minimum allowed value exceeded"),
		MAX_NUMBER_OF_ITERATIONS_EXCEEDED("Maximum number of iterations exceeded");

		private String message;

		private SolverStatus(String message) {
			this.message = message;
		}

		@Override
		public String toString() {
			return message;
		}
	}

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

	public SolverStatus getSolverStatus() {
		return solverStatus;
	}

	void setSolverStatus(SolverStatus solverStatus) {
		this.solverStatus = solverStatus;
	}

	public double getError() {
		return error;
	}

	void setError(double error) {
		this.error = error;
	}
}

