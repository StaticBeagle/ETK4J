package com.wildbitsfoundry.etk4j.math.optimize.solvers;

public class NewtonSolverResults<T> {
	private T value;
	private int numberOfIterations;
	private NewtonSolverStatus solverStatus;
	private double error;

	public enum NewtonSolverStatus {
		CONVERGED,
		MAX_NUMBER_OF_ITERATIONS_EXCEEDED,
		DERIVATIVE_WAS_ZERO,
		TOLERANCE_WAS_REACHED, SECOND_INITIAL_GUESS_EQUAL_TO_FIRST_INITIAL_GUESS
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

	public NewtonSolverStatus getSolverStatus() {
		return solverStatus;
	}

	void setSolverStatus(NewtonSolverStatus newtonSolverStatus) {
		this.solverStatus = newtonSolverStatus;
	}

	public double getError() {
		return error;
	}

	void setError(double error) {
		this.error = error;
	}
}

