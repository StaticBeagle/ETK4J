package examples;

import com.wildbitsfoundry.etk4j.math.optimize.solver.NewtonRaphson;
import com.wildbitsfoundry.etk4j.math.optimize.solver.SolverResults;

public class UnivariateSolverExample {
	public static void main(String[] args) {

		System.out.println("Solve 3 * x^3 - x^2 + 2");
		System.out.println("------------------------------------");
		System.out.printf("Solution with first derivative: ");
		SolverResults<Double> nr = new NewtonRaphson(x ->  x * x * x - x * x + 2, -20)
				.derivative(x -> 3 * x * x - 2 * x)
				.absTolerance(1e-9)
				.relTolerance(0.0)
				.iterationLimit(100)
				.solve();
		System.out.println(nr);
		System.out.println("------------------------------------");
		System.out.printf("Solution with first and second derivative (Halley's method): ");
		SolverResults<Double> nr2 = new NewtonRaphson(x ->  x * x * x - x * x + 2, -20)
				.derivative(x -> 3 * x * x - 2 * x)
				.secondDerivative(x -> 6 * x)
				.solve();
		System.out.println(nr2);
		System.out.println("------------------------------------");
		System.out.printf("Solution with no derivatives (Secant method): ");
		SolverResults<Double> nr3 = new NewtonRaphson(x ->  x * x * x - x * x + 2, -20)
				.solve();
		System.out.println(nr3);
	}
}
