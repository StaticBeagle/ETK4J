package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

/**
 * Schur implements the Schur decomposition of a matrix.  Specifically,
 * given a square matrix A, there is a unitary matrix U such that
 * <pre>
 *      T = U^H AU
 * </pre>
 * is upper triangular.  Schur represents T as a Zutmat and U as a Zmat.
 *
 * @author G. W. Stewart
 * @version Pre-alpha
 */
public class ComplexSchurDecompositionDense {

    /**
     * The upper triangular matrix.
     */
    public ComplexMatrixDense T;

    /**
     * The unitary matrix.
     */
    public ComplexMatrixDense U;

    /**
     * Limits the number of interations in the QR algorithm
     */
    public static int MAXITER = 30;

    /**
     * Creats a Schur decomposition from a square Zmat.
     *
     * @param A The Zmat whose Schur decomposition is to be computed
     *          Thrown for nonsquare matrix.<br>
     *          Thrown for maximum iteration count exceeded.
     */
    public ComplexSchurDecompositionDense(ComplexMatrixDense A) {

        int i, il, iter, iu, k;
        double d, sd, sf;
        Complex b, c, disc, kappa,
                p, q, r, r1, r2, s, z1, z2;
        ComplexPlaneRotationDense P = new ComplexPlaneRotationDense();

        if (A.getRowCount() != A.getColumnCount()) {
//            throw new JampackException // TODO
//                    ("Nonsquare matrix");
        }

        /* Reduce to Hessenberg form and set up T and U */

        ComplexHessembergDecompositionDense H = new ComplexHessembergDecompositionDense(A);
        T = new ComplexMatrixDense(H.H);
        U = H.U;

        iu = T.getRowCount();
        iter = 0;
        while (true) {
            // Locate the range in which to iterate.

            while (iu > 1) {
                d = T.get(iu - 1, iu - 1).norm1() + T.get(iu - 2, iu - 2).norm1();
                sd = T.get(iu - 1, iu - 2).norm1();
                if (sd >= 1.0e-16 * d) break;
                T.set(iu - 1, iu - 2, new Complex());
                iter = 0;
                iu = iu - 1;
            }
            if (iu == 1) break;

            iter = iter + 1;
            if (iter >= MAXITER) {
//                throw new JampackException
//                        ("Maximum number of iterations exceeded."); TODO
            }
            il = iu - 1;
            while (il > 1) {
                d = T.get(il - 1, il - 1).norm1() + T.get(il - 2, il - 2).norm1();
                sd = T.get(il - 1, il - 2).norm1();
                if (sd < 1.0e-16 * d) break;
                il = il - 1;
            }
            if (il != 1) {
                T.set(il - 1, il - 2, new Complex());
            }

            // Compute the shift.

            p = T.get(iu - 2, iu - 2).copy();
            q = T.get(iu - 2, iu - 1).copy();
            r = T.get(iu - 1, iu - 2).copy();
            s = T.get(iu - 1, iu - 1).copy();

            sf = p.norm1() + q.norm1() + r.norm1() + s.norm1();
            p.divideEquals(sf);
            q.divideEquals(sf);
            r.divideEquals(sf);
            s.divideEquals(sf);

            z1 = p.multiply(s);
            z2 = r.multiply(q);
            c = z1.subtract(z2);
            b = p.add(s);

            z1 = b.multiply(b);
            z2 = c.multiply(4);
            disc = z1.subtract(z2).sqrt();

            r1 = b.add(disc).divide(2);
            r2 = b.subtract(disc).divide(2);
            if (r1.norm1() > r2.norm1()) {
                r2 = c.divide(r1);
            } else {
                r1 = c.divide(r2);
            }
            z1 = r1.subtract(s);
            z2 = r2.subtract(s);
            if (z1.norm1() < z2.norm1()) {
                kappa = r1.multiply(sf);
            } else {
                kappa = r2.multiply(sf);
            }
            System.out.println(kappa);
            // Perform the QR step.
            p = T.get(il - 1, il - 1).subtract(kappa).copy();
            q = T.get(il, il - 1).copy();
            ComplexPlaneRotationDense.genc(p, q, P);
            for (i = il; i < iu; i++) {
                ComplexPlaneRotationDense.pa(P, T, i - 1, i, i - 1, T.getColumnCount() - 1);
                ComplexPlaneRotationDense.aph(T, P, 0, Math.min(i + 3, iu) - 1, i - 1, i);
                ComplexPlaneRotationDense.aph(U, P, 0, U.getRowCount() - 1, i - 1, i);
                if (i != iu - 1) {
                    ComplexPlaneRotationDense.genc(T, i, i + 1, i - 1, P);
                }
            }
        }
    }
}
