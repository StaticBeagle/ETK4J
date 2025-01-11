package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;

/**
 * Zhess implements the unitary reduction to Hessenberg form
 * by a unitary similarity transformation. Specifically, given
 * a square matrix A, there is a unitary matrix U such that
 * <pre>
 *      H = U^H AU
 * </pre>
 * is upper Hessenberg.
 * Zhess represents U and H as Zmats.
 *
 * @version Pre-alpha
 * @author G. W. Stewart
 */
public class HessembergDecompositionDense {

    /** The upper Hessenberg matrix */
    public MatrixDense H;

    /** The unitary matrix */
    public MatrixDense U;

    /** Creates a Zhess from a square Zmat. Throws a
     * JampackException for nonsquare matrx.
     *
     * @param A A Zmat
     * Thrown if A is not square.
     */
    public HessembergDecompositionDense(MatrixDense A) {

        if (A.getRowCount() != A.getColumnCount()) {
            //throw new JampackException("Matrix not square"); TODO
        }

        H = new MatrixDense(A);
        U = MatrixDense.Factory.identity(H.getRowCount());

        double[] work = new double[H.getRowCount()];

        for (int k = 0; k < H.getColumnCount() - 2; k++) {
            double[] u = HouseholderTransformations.genc(H, k + 1, H.getRowCount() - 1, k);
            HouseholderTransformations.ua(u, H, k + 1, H.getRowCount() - 1, k + 1, H.getColumnCount() - 1, work);
            HouseholderTransformations.au(H, u, 0, H.getRowCount() - 1, k + 1, H.getColumnCount() - 1, work);
            HouseholderTransformations.au(U, u, 0, U.getRowCount() - 1, k + 1, U.getColumnCount() - 1, work);
        }
    }

    public MatrixDense getH() {
        return H;
    }

    public MatrixDense getU() {
        return U;
    }
}
