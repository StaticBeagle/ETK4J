package com.wildbitsfoundry.etk4j.math.linearalgebra;

/**
 * This class implements the Gauss-Seidel iterative method to solve a system of equations.
 */
public class GaussSeidelSolver extends SuccessiveOverRelaxationSolver {

    public GaussSeidelSolver(MatrixDense A, double[] b) {
        super(A, b, 1.0);
    }
}
