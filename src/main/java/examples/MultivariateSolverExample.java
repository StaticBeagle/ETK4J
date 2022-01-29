package examples;

import com.wildbitsfoundry.etk4j.math.functions.MultivariateFunction;
import com.wildbitsfoundry.etk4j.math.solvers.multivariate.NewtonRaphson;
import com.wildbitsfoundry.etk4j.math.solvers.multivariate.SolverResults;

// TODO refactor multivariate Newton
public class MultivariateSolverExample {

	public static void main(String[] args) {
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

		System.out.printf("Solution with pre-computed Jacobian%n------------------------------------%n");
		SolverResults sol1 = new NewtonRaphson(functions, jacobian, initialguess)
				.setAbsTolerance(1e-9)
				.setRelTolerance(1e-6)
				.setIterationLimit(100)
				.solve();
		System.out.println(sol1);

		System.out.printf("%nSolution with auto-computed Jacobian%n------------------------------------%n");
		SolverResults sol2 = new NewtonRaphson(functions, initialguess)
				.setAbsTolerance(1e-9)
				.setRelTolerance(1e-6)
				.setIterationLimit(100)
				.setDifferentiationStepSize(1e-7)
				.solve();
		System.out.println(sol2);
	}
}
