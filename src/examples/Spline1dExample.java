package examples;

import com.wildbitsfoundry.etk4j.util.NumArrays;
import com.wildbitsfoundry.etk4j.math.interpolation.CubicSpline;

import static com.wildbitsfoundry.etk4j.math.interpolation.CubicSpline.newCubicSpline;


public class Spline1dExample {

	public static void main(String[] args) {
		double[] x = { 1, 2, 3, 4 };
		double[] y = { 5, 6, 7, 8 };
		
		// Creates a cubic spline with not-a-knot conditions
		// Not-a-knot conditions require 4 or more input points
		// This is equivalent to 
		// CubicSpline cs = newNotAKnotSpline(x, y);
		CubicSpline cs = newCubicSpline(x, y);

		// Create an array from 1.0 to 4.0 using a step size of 0.2
		double[] xi = NumArrays.linsteps(1.0, 0.2, 4.0);
		for (double xii : xi) {
			System.out.printf("y(%.4f) = %.4f%n", xii, cs.evaluateAt(xii));
		}

	}
}
