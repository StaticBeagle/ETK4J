package com.wildbitsfoundry.etk4j.math.calculus;

import com.wildbitsfoundry.etk4j.math.functions.MultivariateFunction;

public class CalculusUtil {
    private CalculusUtil() {}

    /**
     * Compute the Jacobian matrix using a 5 point stencil centered difference.
     * @param functions The multi-dimensional function to evaluate.
     * @param x The point at which the Jacobian is computed.
     * @param h Small perturbation for finite differences.
     * @return The Jacobian matrix as a 2D array.
     */
    public static double[][] computeJacobian(MultivariateFunction[] functions, double[] x, double h) {
        int n = x.length; // Number of variables
        int m = functions.length; // Number of equations
        double[][] jacobian = new double[m][n];
        for(int i = 0; i < m; i++) {
            for(int j = 0; j < n; j++) {
                jacobian[i] = PartialDerivatives.centeredDifference5Points(functions[i], x, h);
            }
        }
        return jacobian;
    }
}
