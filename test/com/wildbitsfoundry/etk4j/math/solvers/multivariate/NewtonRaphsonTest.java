package com.wildbitsfoundry.etk4j.math.solvers.multivariate;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.junit.Test;

import com.wildbitsfoundry.etk4j.math.functions.MultivariateFunction;

public class NewtonRaphsonTest {
	
	/*
	 * Solve: 
	 *  f1(x0, x1, x2) = 3 * x0 - cos(x1 * x2) - 1.5 = 0;
	 *  f2(x0, x1, x2) = 4 * x0^2 - 625 * x1^2 + 2 * x2 - 1 = 0;
	 *  f3(x0, x1, x2) = 20 * x2 + exp(-x0 * x1) + 9 = 0;
	 */

	// Define functions
	MultivariateFunction f1 = x -> {
		return 3 * x[0] - Math.cos(x[1] * x[2]) - 1.5;
	};
	MultivariateFunction f2 = x -> {
		return 4 * x[0] * x[0] - 625.0 * x[1] * x[1] + 2.0 * x[2] - 1.0;
	};
	MultivariateFunction f3 = x -> {
		return 20 * x[2] + Math.exp(-x[0] * x[1]) + 9.0;
	};

	// Define partial derivatives
	MultivariateFunction df1x1 = x -> {
		return 3;
	};
	MultivariateFunction df1x2 = x -> {
		return x[2] * Math.sin(x[1] * x[2]);
	};
	MultivariateFunction df1x3 = x -> {
		return x[1] * Math.sin(x[1] * x[2]);
	};

	MultivariateFunction df2x1 = x -> {
		return 8 * x[0];
	};
	MultivariateFunction df2x2 = x -> {
		return -1250.0 * x[1];
	};
	MultivariateFunction df2x3 = x -> {
		return 2.0;
	};

	MultivariateFunction df3x1 = x -> {
		return -x[1] * Math.exp(x[0] * x[1]);
	};

	MultivariateFunction df3x2 = x -> {
		return -x[0] * Math.exp(x[0] * x[1]);
	};

	MultivariateFunction df3x3 = x -> {
		return 20.0;
	};

	// Define Jacobian
	MultivariateFunction[][] jacobian = { { df1x1, df1x2, df1x3 }, { df2x1, df2x2, df2x3 },
			{ df3x1, df3x2, df3x3 } };

	MultivariateFunction[] functions = { f1, f2, f3 };
	double[] initialguess = { 1.0, 1.0, 1.0 };

	
	@Test
	public void testPreComputedJacobian() {
		SolverResults sol = new NewtonRaphson(functions, jacobian, initialguess).setAbsTolerance(1e-9).setRelTolerance(1e-6)
				.setIterationLimit(100).solve();
		assertArrayEquals(new double[] { 0.8332816138167558, 0.03533461613948315, -0.49854927781116914}, sol.Solution, 1e-12);
	}
	
	@Test
	public void testSingularJacobian() {
		MultivariateFunction df1x1 = x -> { return 0.0; };
		MultivariateFunction df1x2 = x -> { return 0.0; };
		MultivariateFunction df1x3 = x -> { return 0.0; };
		MultivariateFunction[][] jacobian = { { df1x1, df1x2, df1x3 }, { df2x1, df2x2, df2x3 },
				{ df3x1, df3x2, df3x3 } };

		MultivariateFunction[] functions = { f1, f2, f3 };
		double[] initialguess = { 1.0, 1.0, 1.0 };

		SolverResults sol = new NewtonRaphson(functions, jacobian, initialguess).setAbsTolerance(1e-9).setRelTolerance(1e-6)
				.setIterationLimit(100).solve();
		assertEquals(SolverResults.SolverStatus.JACOBIAN_IS_SINGULAR, sol.Status);
	}
	
	@Test
	public void testMaxValueExceeded() {
		SolverResults sol = new NewtonRaphson(functions, jacobian, initialguess).setAbsTolerance(1e-9).setRelTolerance(1e-6)
				.setIterationLimit(100).setMaxAbsAllowedValue(0.0).solve();
		assertEquals(SolverResults.SolverStatus.MAX_ABS_VALUE_EXCEEDED, sol.Status);
	}
	
	@Test
	public void testMaxIterationsExceeded() {
		SolverResults sol = new NewtonRaphson(functions, jacobian, initialguess).setAbsTolerance(1e-9).setRelTolerance(1e-6)
				.setIterationLimit(1).solve();
		assertEquals(SolverResults.SolverStatus.ITERATION_LIMIT_EXCEEDED, sol.Status);
	}
	
	@Test
	public void testMaxAbsErrorScheme() {
		SolverResults sol = new NewtonRaphson(functions, jacobian, initialguess).setAbsTolerance(1e-9).setRelTolerance(1e-6)
				.setIterationLimit(100).setErrorEstimationScheme(ErrorEstimationScheme.MAX_ABS_ERROR).solve();
		assertArrayEquals(new double[] { 0.8332816138167558, 0.03533461613948315, -0.49854927781116914}, sol.Solution, 1e-12);
	}
	
	@Test
	public void testSumAbsErrorScheme() {
		SolverResults sol = new NewtonRaphson(functions, jacobian, initialguess).setAbsTolerance(1e-9).setRelTolerance(1e-6)
				.setIterationLimit(100).setErrorEstimationScheme(ErrorEstimationScheme.SUM_ABS_ERROR).solve();
		assertArrayEquals(new double[] { 0.8332816138167558, 0.03533461613948315, -0.49854927781116914}, sol.Solution, 1e-12);
	}
	
	@Test
	public void testCalculatedJacobian() {
		SolverResults sol = new NewtonRaphson(functions, initialguess).setAbsTolerance(1e-9).setRelTolerance(1e-6)
				.setIterationLimit(100).setDifferentiationStepSize(1e-7).solve();
		assertArrayEquals(new double[] { 0.8332816138167558, 0.03533461613948315, -0.49854927781116914}, sol.Solution, 1e-12);
	}
}
