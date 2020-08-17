package examples;

import static com.wildbitsfoundry.etk4j.math.interpolation2d.Spline2d.newBilinearSpline;

import com.wildbitsfoundry.etk4j.math.interpolation2d.Spline2d;

public class Spline2DifferentDimensionsExample {
	
	public static void main(String[] args) {
		createPlane();
	}
	


	public static void createPlane() {
		double[] x = { 0, 1, 2, 3, 4, 5, 6, 7 };
		double[] y = { 1, 2 };

		double[][] z = { { 1, 2, 3, 4, 5, 6, 7, 8 }, { 2, 4, 6, 8, 10, 12, 14, 16 } };

		System.out.printf("%n:: Test Bilinear spline%n");
		Spline2d sp = newBilinearSpline(x, y, z);
		Spline2dExample.printSpline2d(sp, x, y);
	}
}

