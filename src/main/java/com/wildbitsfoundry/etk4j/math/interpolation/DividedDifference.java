package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;

public class DividedDifference implements UnivariateFunction {

    private final double[][] dividedDifferenceTable;
    private final double[] x;

    public DividedDifference(double[] x, double[] y) {
        dividedDifferenceTable = computeDividedDifferenceTable(x, y);
        this.x = x;
    }

    @Override
    public double evaluateAt(double x) {
        return newtonPolynomial(this.x, dividedDifferenceTable, x);
    }

    // Method to compute the divided difference table
    public static double[][] computeDividedDifferenceTable(double[] x, double[] y) {
        int n = x.length;
        double[][] table = new double[n][n];

        // Initialize the first column with y-values
        for (int i = 0; i < n; i++) {
            table[i][0] = y[i];
        }

        // Compute the divided differences iteratively
        for (int j = 1; j < n; j++) {
            for (int i = 0; i < n - j; i++) {
                table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (x[i + j] - x[i]);
            }
        }

        return table;
    }

    // Method to construct the Newton interpolating polynomial
    public static double newtonPolynomial(double[] x, double[][] table, double value) {
        int n = x.length;
        double result = table[0][0];
        double product;

        for (int i = 1; i < n; i++) {
            product = 1.0;
            for (int j = 0; j < i; j++) {
                product *= (value - x[j]);
            }
            result += table[0][i] * product;
        }

        return result;
    }
}
