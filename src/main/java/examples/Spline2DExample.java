package examples;

import com.wildbitsfoundry.etk4j.math.functions.BivariateFunction;
import com.wildbitsfoundry.etk4j.math.interpolation.BiCubicSpline;
import com.wildbitsfoundry.etk4j.math.interpolation.BiLinearSpline;
import com.wildbitsfoundry.etk4j.math.interpolation.Spline2D;

public class Spline2DExample {

	public static void main(String[] args) {
		createSurface();
	}
	
	public static void createSurface() {
		double[] x = { 1, 2, 3, 4, 5, 6, 7, 8 };
		double[] y = { 1, 2, 3, 4, 5, 6, 7, 8 };

		double[][] z = { { 1, 4, 9, 16, 25, 36, 49, 64 }, { 4, 16, 36, 64, 100, 144, 196, 256 },
				{ 9, 36, 81, 144, 225, 324, 441, 576 }, { 16, 64, 144, 256, 400, 576, 784, 1024 },
				{ 25, 100, 225, 400, 625, 900, 1225, 1600 }, { 36, 144, 324, 576, 900, 1296, 1764, 2304 },
				{ 49, 196, 441, 784, 1225, 1764, 2401, 3136 }, { 64, 256, 576, 1024, 1600, 2304, 3136, 4096 } };

		System.out.printf(":: Function z(x, y) = x^2 * y^2%n%n");

		System.out.println(":: Test Bicubic spline");
		Spline2D sp = BiCubicSpline.newBicubicSpline(x, y, z);
		printSpline2D(sp, x, y);

		System.out.printf("%n:: Test Bilinear spline%n");
		sp = BiLinearSpline.newBilinearSpline(x, y, z);
		printSpline2D(sp, x, y);
	}

	public static void printSpline2D(BivariateFunction sp, double[] x, double[] y) {

		final int rows = y.length;
		final int cols = x.length;
		System.out.print("   x |");
		for (int i = 0; i < cols; ++i) {
			System.out.printf("%7.2f ", x[i]);
		}
		System.out.printf("%ny    |                           z(x, y)%n");
		System.out.printf("-----|---------------------------------------------------------------%n");
		for (int i = 0; i < rows; ++i) {
			System.out.printf("%.2f |", y[i]);
			for (int j = 0; j < cols; ++j) {
				System.out.printf("%7.2f ", sp.evaluateAt(x[j], y[i]));
			}
			System.out.println();
		}
		System.out.println();

		for (int i = 1; i < rows; ++i) {
			double pi = i + 0.5;
			System.out.printf("z(%.1f, %.1f) = %.4f%n", pi, pi, sp.evaluateAt(pi, pi));
		}
	}
}
