package com.wildbitsfoundry.etk4j.math.linearalgebra;

// TODO javadocs
public final class LeastSquaresSolver {

    private LeastSquaresSolver() {}

    public static MatrixDense solve(MatrixDense A, double[] b) {
        return new QRDecompositionDense(A).solve(new MatrixDense(b, b.length));
    }

    public static MatrixDense solve(MatrixDense A, MatrixDense b) {
        return new QRDecompositionDense(A).solve(b);
    }
}
