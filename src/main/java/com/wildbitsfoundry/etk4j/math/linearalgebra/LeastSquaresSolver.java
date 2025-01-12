package com.wildbitsfoundry.etk4j.math.linearalgebra;

/**
 * This class solves an over determined (more solutions than unknowns) in a least square sense
 */
public final class LeastSquaresSolver {

    private LeastSquaresSolver() {}

    /**
     * Solves the linear system A*x = b using the {@link QRDecompositionDense} algorithm
     * @param A The coefficient Matrix
     * @param b The solution vector
     * @return x = A<sup>-1</sup>*b
     */
    public static MatrixDense solve(MatrixDense A, double[] b) {
        return new QRDecompositionDense(A).solve(new MatrixDense(b, b.length));
    }

    /**
     * Solves the linear system A*x = b using the {@link QRDecompositionDense} algorithm
     * @param A The coefficient Matrix
     * @param b The solution Matrix
     * @return x = A<sup>-1</sup>*b
     */
    public static MatrixDense solve(MatrixDense A, MatrixDense b) {
        return new QRDecompositionDense(A).solve(b);
    }
}
