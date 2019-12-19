package com.wildbitsfoundry.etk4j.statistics.regression;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;
import com.wildbitsfoundry.etk4j.util.NumArrays;

class MultiplePolynomialRegression extends MultivariateRegression {
	public MultiplePolynomialRegression(double[][] x, double[] y, int[] n) {

		// int order = (m + 1) * (n + 1); // ... * (nx + 1)
		int[] degrees = n;
		int dims = x.length;
		// double[][] vars = new double[dims][];
		// vars[0] = x;
		// vars[1] = y;
		// ... vars[n]

		int k = dims - 1;
		double[][] vandermonde = new double[y.length][];
		for (int i = 0; i < y.length; ++i) {
			vandermonde[i] = buildPowerSeriesDescending(x[k][i], degrees[k]);
			for (int j = k - 1; j >= 0; --j) {
				double[] tmp = buildPowerSeriesDescending(x[j][i], degrees[j]);
				vandermonde[i] = NumArrays.kron(vandermonde[i], tmp);
			}
		}
		Matrix X = new Matrix(vandermonde);
		this.doRegression(X, x, y);
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

	public static void main(String[] args) {
		// double[][] x = { { 10, 20, 40, 80, 160, 200 }, { 20, 40, 15, 100, 23,
		// 18 }, { 15, 40, 62, 76, 44, 18 } };
		// double[] y = { 243, 483, 508, 1503, 1764, 2129 };
//		double[][] x = { { 2, 5, 8, 15, 30, 45, 50, 60, 70, 80 }, { 15, 229, 3480, 15, 229, 3480, 5, 15, 123, 150 },
//				{ 3, 4, 5, 8, 9, 10, 11, 55, 80, 160 } };
		double[][] x = { { 2, 5, 8, 15, 30, 45, 50, 60, 70, 80 }, { 15, 229, 3480, 15, 229, 3480, 5, 15, 123, 150 }};
		double[] y = { 0.5, 105, 1776, 6, 110, 4300, 10000, 12000, 15000, 20000 };
		MultiplePolynomialRegression reg = new MultiplePolynomialRegression(x, y, new int[] { 2, 2 });
		System.out.println(Arrays.toString(reg.beta()));

	}
}
