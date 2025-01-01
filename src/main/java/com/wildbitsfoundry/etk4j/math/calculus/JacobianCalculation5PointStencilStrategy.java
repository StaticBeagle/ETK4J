package com.wildbitsfoundry.etk4j.math.calculus;

import com.wildbitsfoundry.etk4j.math.functions.MultivariateFunction;

public class JacobianCalculation5PointStencilStrategy implements JacobianCalculationStrategy {

    @Override
    public double[][] calculateJacobian(MultivariateFunction[] functions, double[] x, Object... params) {
        int n = x.length; // Number of variables
        int m = functions.length; // Number of equations
        double[][] jacobian = new double[m][n];
        for(int i = 0; i < m; i++) {
            for(int j = 0; j < n; j++) {
                jacobian[i] = PartialDerivatives.centeredDifference5Points(functions[i], x, (double) params[0]);
            }
        }
        return jacobian;
    }
}
