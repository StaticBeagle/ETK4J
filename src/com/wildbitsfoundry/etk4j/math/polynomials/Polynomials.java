package com.wildbitsfoundry.etk4j.math.polynomials;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrices;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;
import com.wildbitsfoundry.etk4j.util.ArrayUtils;

public final class Polynomials {
	private Polynomials() {
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
