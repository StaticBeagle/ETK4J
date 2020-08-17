package examples;

import com.wildbitsfoundry.etk4j.math.solvers.univariate.NewtonRaphson;
import com.wildbitsfoundry.etk4j.math.solvers.univariate.SolverResults;

public class UnivariateSolverExample {
	public static void main(String[] args) {
		/*
		 * Find the square root of 2
		 * Solve: 
		 *  f(x) = 2 - x^2 = 0;
		 */

		/* @formatter:off */
		System.out.printf("Solution with pre-computed Jacobian%n------------------------------------%n");
		SolverResults sr1 = new NewtonRaphson(x -> 2 - x * x, x -> -2 * x, 1)
				.absTolerance(1e-9)
				.relTolerance(1e-6)
				.iterationLimit(100)
				.solve();
		System.out.println(sr1);

		System.out.printf("%nSolution with auto-computed Jacobian%n------------------------------------%n");
		SolverResults sr2 = new NewtonRaphson(x -> 2 - x * x, 1)
				.absTolerance(1e-9)
				.relTolerance(1e-6)
				.iterationLimit(100)
				.differentiationStepSize(1e-7)
				.solve();
		System.out.println(sr2);
	}
}
