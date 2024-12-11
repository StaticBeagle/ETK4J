package examples;

import static com.wildbitsfoundry.etk4j.math.interpolation.Spline2D.newBilinearSpline;

import com.wildbitsfoundry.etk4j.math.interpolation.Spline2D;

public class Spline2DifferentDimensionsExample {
	
	public static void main(String[] args) {
		createPlane();
	}

	public static void createPlane() {
		double[] x = { 0, 1, 2, 3, 4, 5, 6, 7 };
		double[] y = { 1, 2, 3 };

		double[][] z = { { 1, 2, 3, 4, 5, 6, 7, 8 }, { 2, 4, 6, 8, 10, 12, 14, 16 }, {3, 6, 9, 12, 15, 18, 21, 24} };

		System.out.printf("%n:: Test Bilinear spline%n");
		Spline2D sp = newBilinearSpline(x, y, z);
		Spline2DExample.printSpline2D(sp, x, y);
	}
}

