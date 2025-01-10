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
public class ComplexHessembergDecompositionDense {

    /** The upper Hessenberg matrix */
    public ComplexMatrixDense H;

    /** The unitary matrix */
    public ComplexMatrixDense U;

    /** Creates a Zhess from a square Zmat. Throws a
     * JampackException for nonsquare matrx.
     *
     * @param A A Zmat
     * Thrown if A is not square.
     */
    public ComplexHessembergDecompositionDense(ComplexMatrixDense A) {

        if (A.getRowCount() != A.getColumnCount()) {
            //throw new JampackException("Matrix not square"); TODO
        }

        H = new ComplexMatrixDense(A);
        U = ComplexMatrixDense.Factory.identity(H.getRowCount());

        Complex[] work = ComplexArrays.zeros(H.getRowCount());

        for (int k = 0; k < H.getColumnCount() - 2; k++) {
            Complex[] u = ComplexHouseholderTransformations.genc(H, k + 1, H.getRowCount() - 1, k);
            ComplexHouseholderTransformations.ua(u, H, k + 1, H.getRowCount() - 1, k + 1, H.getColumnCount() - 1, work);
            ComplexHouseholderTransformations.au(H, u, 0, H.getRowCount() - 1, k + 1, H.getColumnCount() - 1, work);
            ComplexHouseholderTransformations.au(U, u, 0, U.getRowCount() - 1, k + 1, U.getColumnCount() - 1, work);
        }
    }
}
