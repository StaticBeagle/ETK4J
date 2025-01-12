package com.wildbitsfoundry.etk4j.math.linearalgebra;

/**
 * HessembergDecompositionDense implements the unitary reduction to Hessenberg form
 * by a unitary similarity transformation. Specifically, given
 * a square matrix A, there is a unitary matrix U such that
 * <pre>
 *      H = U^H AU
 * </pre>
 * is upper Hessenberg.
 * HessembergDecompositionDense represents U and H as {@link ComplexMatrix}.
 *
 * @version Pre-alpha
 * @author G. W. Stewart
 */
public class HessembergDecompositionDense {

    /** The upper Hessenberg matrix */
    public MatrixDense H;

    /** The unitary matrix */
    public MatrixDense U;

    /** Creates a Hessemberg Matrix  from a square Complex Matrix. Throws a
     * {@link NonSquareMatrixException} for nonsquare matrx.
     *
     * @param A A {@link ComplexMatrix}
     * @throws NonSquareMatrixException Thrown if A is not square.
     */
    public HessembergDecompositionDense(MatrixDense A) {

        if (A.getRowCount() != A.getColumnCount()) {
            throw new NonSquareMatrixException("Matrix not square");
        }

        H = new MatrixDense(A);
        U = MatrixDense.Factory.identity(H.getRowCount());

        double[] work = new double[H.getRowCount()];

        for (int k = 0; k < H.getColumnCount() - 2; k++) {
            double[] u = HouseholderTransformationsDense.genc(H, k + 1, H.getRowCount() - 1, k);
            HouseholderTransformationsDense.ua(u, H, k + 1, H.getRowCount() - 1, k + 1, H.getColumnCount() - 1, work);
            HouseholderTransformationsDense.au(H, u, 0, H.getRowCount() - 1, k + 1, H.getColumnCount() - 1, work);
            HouseholderTransformationsDense.au(U, u, 0, U.getRowCount() - 1, k + 1, U.getColumnCount() - 1, work);
        }
    }

    /**
     * The upper Hessemberg Matrix
     * @return The upper Hessemberg Matrix
     */
    public MatrixDense getH() {
        return H;
    }

    /**
     * THe unitary Matrix
     * @return The unitary Matrix
     */
    public MatrixDense getU() {
        return U;
    }
}
