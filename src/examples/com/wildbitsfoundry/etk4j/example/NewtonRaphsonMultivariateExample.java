package examples.com.wildbitsfoundry.etk4j.example;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

import com.wildbitsfoundry.etk4j.math.functions.MultivariateFunction;
import com.wildbitsfoundry.etk4j.math.solvers.multivariate.NewtonRaphson;

public class NewtonRaphsonMultivariateExample {
	
	public static void main(String[] args) {
		/*
		 * Solve: 
		 *   f1(x0, x1, x2) = 3 * x0 - cos(x1 * x2) - 1.5 = 0
		 *   f2(x0, x1, x2) = 4 * x0^2 - 625 * x1^2 + 2 * x2 - 1 = 0
		 *   f3(x0, x1, x2) = 20 * x2 + exp(-x0 * x1) + 9 = 0
		 */
		// Uses lambdas
		System.out.println(":: Using functional interface:");
		usingFunctionalInterface();
		
		// Uses Anonymous instantiations of MultivariateFunction
		// In general, this implementation might perform faster 
		// than its functional counterpart
		System.out.printf("%n:: Using multivariate interface:%n");
		usingMultivariateInterface();
	}
	
	public static void usingFunctionalInterface() {
		double[] initialguess = { 1d, 2d, 3d };

		List<Function<Double[], Double>> functions = new ArrayList<>();
		functions.add(x -> 3 * x[0] - Math.cos(x[1] * x[2]) - 1.5);
		functions.add(x -> 4 * x[0] * x[0] - 625.0 * x[1] * x[1] + 2.0 * x[2] - 1.0);
		functions.add(x -> 20 * x[2] + Math.exp(-x[0] * x[1]) + 9.0);

		List<Function<Double[], Double[]>> jacobian = new ArrayList<>();
		jacobian.add(x -> new Double[] { 3.0, x[2] * Math.sin(x[1] * x[2]), x[1] * Math.sin(x[1] * x[2]) });
		jacobian.add(x -> new Double[] { 8 * x[0], -1250.0 * x[1], 2.0 });
		jacobian.add(x -> new Double[] { -x[1] * Math.exp(x[0] * x[1]), -x[0] * Math.exp(x[0] * x[1]), 20.0 });

		System.out.printf("Solution with pre-computed Jacobian%n------------------------------------%n");
		double[] sol1 = new NewtonRaphson(functions, jacobian, initialguess)
				.absTolerance(1e-9)
				.relTolerance(1e-6)
				.iterationLimit(100)
				.solve();
		for (double val : sol1) {
			System.out.printf("%.6g%n", val);
		}

		System.out.printf("%nSolution with auto-computed Jacobian%n------------------------------------%n");
		double[] sol2 = new NewtonRaphson(functions, initialguess)
				.absTolerance(1e-9)
				.relTolerance(1e-6)
				.iterationLimit(100)
				.differentiationStepSize(0.001)
				.solve();
		for (double val : sol2) {
			System.out.printf("%.6g%n", val);
		}
	}

	public static void usingMultivariateInterface() {
		MultivariateFunction f1 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return 3 * x[0] - Math.cos(x[1] * x[2]) - 1.5;
			}
		};

		MultivariateFunction f2 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return 4 * x[0] * x[0] - 625.0 * x[1] * x[1] + 2.0 * x[2] - 1.0;
			}
		};

		MultivariateFunction f3 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return 20 * x[2] + Math.exp(-x[0] * x[1]) + 9.0;
			}
		};

		// Jacobian
		MultivariateFunction df1x1 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return 3;
			}
		};

		MultivariateFunction df1x2 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return x[2] * Math.sin(x[1] * x[2]);
			}
		};

		MultivariateFunction df1x3 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return x[1] * Math.sin(x[1] * x[2]);
			}
		};

		MultivariateFunction df2x1 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return 8 * x[0];
			}
		};

		MultivariateFunction df2x2 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return -1250.0 * x[1];
			}
		};

		MultivariateFunction df2x3 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return 2.0;
			}
		};

		MultivariateFunction df3x1 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return -x[1] * Math.exp(x[0] * x[1]);
			}
		};

		MultivariateFunction df3x2 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return -x[0] * Math.exp(x[0] * x[1]);
			}
		};

		MultivariateFunction df3x3 = new MultivariateFunction() {

			@Override
			public double evaluateAt(double... x) {
				return 20.0;
			}
		};

		MultivariateFunction[][] jacobian = new MultivariateFunction[][] { { df1x1, df1x2, df1x3 },
				{ df2x1, df2x2, df2x3 }, { df3x1, df3x2, df3x3 } };

		MultivariateFunction[] functions = new MultivariateFunction[] { f1, f2, f3 };
		double[] initialguess = new double[] { 1d, 2d, 3d };

		System.out.printf("Solution with pre-computed Jacobian%n------------------------------------%n");
		double[] sol1 = new NewtonRaphson(functions, jacobian, initialguess)
				.absTolerance(1e-9)
				.relTolerance(1e-8)
				.iterationLimit(100)
				.solve();
		for (double val : sol1) {
			System.out.printf("%.6g%n", val);
		}

		System.out.printf("%nSolution with auto-computed Jacobian%n------------------------------------%n");
		double[] sol2 = new NewtonRaphson(functions, initialguess)
				.absTolerance(1e-9)
				.relTolerance(1e-8)
				.iterationLimit(100)
				.differentiationStepSize(0.001)
				.solve();
		for (double val : sol2) {
			System.out.printf("%.6g%n", val);
		}
	}
}
