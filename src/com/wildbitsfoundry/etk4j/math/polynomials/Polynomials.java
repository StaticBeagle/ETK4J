package com.wildbitsfoundry.etk4j.math.polynomials;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrices;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;
import com.wildbitsfoundry.etk4j.util.ArrayUtils;

public final class Polynomials {
	private Polynomials() {
	}

	public static Complex[] quadraticFormula(double a, double b, double c) {
		double k = 1 / (2 * a);
		double dis = b * b - 4 * a * c;
		Complex r1;
		Complex r2;
		if (dis < 0) {
			r1 = new Complex(-b * k, Math.sqrt(-dis) * k);
			r2 = r1.conj();
		} else {
			r1 = new Complex(k * (-b + Math.sqrt(dis)), 0.0);
			r2 = new Complex(k * (-b - Math.sqrt(dis)), 0.0);
		}
		return new Complex[] { r1, r2 };
	}

	public static Complex[] cubicFormula(double a, double b, double c, double d) {

		double a2 = a * a;
		double b2 = b * b;
		double c2 = c * c;
		double abc = a * b * c;
		double disc = 18.0 * abc * d - 4.0 * b2 * b * d + b2 * c2 - 4.0 * a * c2 * c - 27.0 * a2 * d * d;
		double disc0 = b2 - 3.0 * a * c;
		double disc1 = 2.0 * b2 * b - 9.0 * abc + 27.0 * a2 * d;

		Complex k1 = new Complex(-0.5, 0.866025403784439);
		Complex k2 = new Complex(-0.5, -0.866025403784439);
		Complex r = new Complex(-27.0 * a2 * disc, 0.0).sqrt();
		Complex C = r.add(disc1).divide(2d).pow(1.0 / 3.0);
		double k = -1.0 / (3.0 * a);
		return new Complex[] { C.add(b).add(C.pow(-1).multiply(disc0)).multiply(k),
				C.multiply(k1).add(b).add(C.multiply(k1).pow(-1).multiply(disc0)).multiply(k),
				C.multiply(k2).add(b).add(C.multiply(k2).pow(-1).multiply(disc0)).multiply(k) };
	}
	
	/***
	 * Polynomial fit
	 * <pre>
	 * Finds a polynomial P(x) = c0 + c1*x + * c2*x^2 + ... + cn*x^n 
	 * of degree n that fits the data in y best in a least-square sense.
	 * </pre>
	 * @param x
	 *            Array of points for the independent variable x
	 * @param y
	 *            Array of solutions to y(x)
	 * @param n
	 *            Order of the polynomial
	 * @return Returns a polynomial of degree n fits the data y best in a
	 *         least-square sense
	 */
	public static Polynomial polyfit(double[] x, double[] y, int n) {
		int dim = x.length;
		if (dim != y.length) {
			throw new IllegalArgumentException("x and y dimensions must match!");
		}
		// Building the coefficient matrix
		Matrix A = Matrices.Vandermonde(x, dim, n);
		// Building the solution vector
		Matrix b = new Matrix(y, dim);
		Matrix c = A.solve(b);

		double[] coeffs = new double[n + 1];
		for (int i = 0; i <= n; i++) {
			coeffs[i] = c.get(n - i, 0);
		}
		Polynomial result = new Polynomial(coeffs);
		return result;
	}

	public static void main(String[] args) {
		double[] x = new double[] { 1, 2, 3, 4 };
		double[] y = new double[] { 3, 4, 5, 6 };

		double[] z = new double[] { 5, 6, 7, 8 };

		System.out.println(Arrays.toString(polyfit2D(x, y, z, 2, 2)));
	}
	
	public static double[] polyfit2D(double[] x, double[] y, double[] z, int nx, int ny) {
		double[][] vars = new double[2][];
		vars[0] = x;
		vars[1] = y;
		
		int[] n = new int[] { nx, ny };
		return polyfitN(vars, z, n);
	}
	
	private static double[] polyfitN(double[][] vars, double[] z, int[] n) {

		//int order = (m + 1) * (n + 1); // ... * (nx + 1)
		int[] degrees = new int[] { 1, 1 };
		int dims = vars.length;
//		double[][] vars = new double[dims][];
//		vars[0] = x;
//		vars[1] = y;
		// ... vars[n]

		int k = dims - 1;
		double[][] vandermonde = new double[z.length][];
		for (int i = 0; i < z.length; ++i) {
			vandermonde[i] = buildPowerSeriesDescending(vars[k][i], degrees[k]);
			for (int j = k - 1; j >= 0; --j) {
				double[] tmp = buildPowerSeriesDescending(vars[j][i], degrees[j]);
				vandermonde[i] = ArrayUtils.kron(vandermonde[i], tmp);
			}
		}
		
		Matrix A = new Matrix(vandermonde);
		Matrix b = new Matrix(z, z.length);
		
		Matrix c = A.solve(b);

		return c.getArrayCopy();
	}
	
	private static double[] buildPowerSeriesDescending(double x, int power) {
		double[] series = new double[power + 1];
		int i = 0;
		while (power > 0) {
			series[i++] = Math.pow(x, power--);
		}
		series[i] = 1.0;
		return series;
	}
}
