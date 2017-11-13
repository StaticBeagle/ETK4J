package com.wildbitsfoundry.etk4j.statistics.regression;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;

public class MultipleLinearRegression extends MultivariateRegression {

	public MultipleLinearRegression(double[][] x, double[] y) {

		final int rows = x.length;
		final int cols = x[0].length + 1;
		double[][] xp = new double[rows][cols];
		for (int i = 0; i < rows; ++i) {
			xp[i][0] = 1.0;
			System.arraycopy(x[i], 0, xp[i], 1, cols - 1);
		}
		Matrix X = new Matrix(xp);
		this.doRegression(X, x, y);
	}

	public static void main(String[] args) {
		double[][] x = { { 10, 20 }, { 20, 40 }, { 40, 15 }, { 80, 100 }, { 160, 23 }, { 200, 18 } };
		double[] y = { 243, 483, 508, 1503, 1764, 2129 };
		MultipleLinearRegression reg = new MultipleLinearRegression(x, y);
		System.out.println(Arrays.toString(reg.beta()));


	}
}
