package com.wildbitsfoundry.etk4j.math.solvers.univariate;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import com.wildbitsfoundry.etk4j.math.solvers.univariate.SolverResults.SolverStatus;

public class NewtonRaphsonTest {
	
	@Test
	public void testPreComputedDerivative() {
		SolverResults sol = new NewtonRaphson(x -> 2 - x * x, x -> -2 * x, 1)
				.absTolerance(1e-9)
				.relTolerance(1e-6)
				.iterationLimit(100)
				.solve();
		assertEquals(1.4142135623730951, sol.Solution, 1e-12);
	}
	
	@Test
	public void testMaxValueExceeded() {
		SolverResults sol = new NewtonRaphson(x -> 2 - x * x, x -> -2 * x, 1)
				.absTolerance(1e-9)
				.relTolerance(1e-6)
				.iterationLimit(100)
				.maxAbsAllowedValue(0.0)
				.solve();
		assertEquals(SolverResults.SolverStatus.MAX_ABS_VALUE_EXCEEDED, sol.Status);
	}
	
	@Test
	public void testMaxIterationsExceeded() {
		SolverResults sol = new NewtonRaphson(x -> 2 - x * x, x -> -2 * x, 1)
				.absTolerance(1e-9)
				.relTolerance(1e-6)
				.iterationLimit(100)
				.iterationLimit(1)
				.solve();
		assertEquals(SolverResults.SolverStatus.ITERATION_LIMIT_EXCEEDED, sol.Status);
	}
	
	@Test
	public void testCalculatedDerivative() {
		SolverResults sol = new NewtonRaphson(x -> 2 - x * x, 1)
				.absTolerance(1e-9)
				.relTolerance(1e-6)
				.iterationLimit(100)
				.differentiationStepSize(1e-7)
				.solve();
		assertEquals(1.4142135623730951, sol.Solution, 1e-12);
	}
	
	@Test
	public void testToString() {
		SolverResults sol = new NewtonRaphson(x -> 2 - x * x, x -> -2 * x, 1)
				.absTolerance(1e-9)
				.relTolerance(1e-6)
				.iterationLimit(100)
				.solve();
		String expected = String.format("Status: SUCCESS%n" +
				"Iterations: 5%n" +
				"Converged: true%n" +
				"Estimated Error: 1.5947243525715749E-12%n" +
				"Solution: 1.41421");
		assertEquals(expected, sol.toString());
	}
	
	@Test
	public void testSolverStatusEnum() {
		SolverStatus status = SolverStatus.valueOf("SUCCESS");
		assertEquals(SolverStatus.SUCCESS, status);
	}
}
