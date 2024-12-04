package com.wildbitsfoundry.etk4j.math.linearalgebra;
// TODO javadocs
public class GaussSeidelSolver extends SuccessiveOverRelaxationSolver {

    public GaussSeidelSolver(MatrixDense A, double[] b) {
        super(A, b, 1.0);
    }
}
