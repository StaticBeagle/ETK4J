package com.wildbitsfoundry.etk4j.math.linearalgebra;

public class GaussSeidel extends SuccessiveOverRelaxation {

    public GaussSeidel(MatrixDense A, double[] b) {
        super(A, b, 1.0);
    }
}
