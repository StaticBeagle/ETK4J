package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.constant.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.signal.filters.MaximumNumberOfIterationsExceededException;

/*
 * Zsvd implements the singular value decomposion of a ComplexMatrixDense.
 * Specifically if X is an mxn matrix with m&gt;=n there are unitary
 * matrices U and V such that
 * <pre>
 *     U^H*X*V = | S |
 *               | 0 |
 * </pre>
 * where S = diag(s1,...,sm) with
 * <pre>
 *     s1 >= s2 >= ... >= sn >=0.
 * </pre>
 * If m&lt;n the decomposition has the form
 * <pre>
 *     U^H*X*V = | S  0 |,
 * </pre>
 * <p>
 * where S is diagonal of order m.  The diagonals of S are the
 * singular values of A.  The columns of U are the left singular
 * vectors of A and the columns of V are the right singular vectors.
 *
 * @author G. W. Stewart
 * @version Pre-alpha
 */
public class ComplexSingularValueDecompositionDense {

    /**
     * Limits the number of iterations in the SVD algorithm
     */
    public static int MAXITER = 30;

    /**
     * The matrix of left singular vectors
     */
    public ComplexMatrixDense U;

    /**
     * The matrix of right singular vectors
     */
    public ComplexMatrixDense V;

    /**
     * The diagonal matrix of singular values
     */

    public MatrixDense S;

    private final int rows;
    private final int cols;
    private final double[] d;
    /*
     Computes the SVD of a ComplexMatrixDense XX.  Throws a JampackException
     if the maximum number of iterations is exceeded.

     @param     XX A ComplexMatrixDense
     @return    The Zsvd of XX
     @exception JampackException
     Thrown if maximimum number of iterations is
     exceeded.<br>
     Passed from below.
     */

    public ComplexSingularValueDecompositionDense(ComplexMatrixDense XX) {

        rows = XX.getRowCount();
        cols = XX.getColumnCount();

        int i, il, iu, iter, j, k, kk, m, mc;

        double as, at, au, axkk, axkk1, dmax, dmin, ds, ea,
                es, shift, ss, t, tre;

        Complex xkk, xkk1, xk1k1, ukj, vik1;

        ComplexPlaneRotationDense P = new ComplexPlaneRotationDense();

        /*    Initialization */

        Complex scale;
        Complex zr;

        ComplexMatrixDense X = new ComplexMatrixDense(XX);

        Complex[] h;
        Complex[] temp = new Complex[Math.max(X.getRowCount(), X.getColumnCount())];

        mc = Math.min(X.getRowCount(), X.getColumnCount());
        d = new double[mc];
        double[] e = new double[mc];

        S = new MatrixDense(mc, mc);
        U = ComplexMatrixDense.Factory.identity(X.getRowCount());
        V = ComplexMatrixDense.Factory.identity(X.getColumnCount());

        m = Math.min(X.getRowCount() - 1, X.getColumnCount() - 1);
      

/*
      Reduction to Bidiagonal form.
*/

        for (k = 0; k <= m; k++) {

            h = ComplexHouseholderTransformationsDense.genc(X, k, X.getRowCount() - 1, k);
            ComplexHouseholderTransformationsDense.ua(h, X, k, X.getRowCount() - 1, k + 1, X.getColumnCount() - 1, temp);
            ComplexHouseholderTransformationsDense.au(U, h, 0, U.getRowCount() - 1, k, U.getColumnCount() - 1, temp);

            if (k != X.getColumnCount() - 1) {
                h = ComplexHouseholderTransformationsDense.genr(X, k, k + 1, X.getColumnCount() - 1);
                ComplexHouseholderTransformationsDense.au(X, h, k + 1, X.getRowCount() - 1, k + 1, X.getColumnCount() - 1, temp);
                ComplexHouseholderTransformationsDense.au(V, h, 0, V.getRowCount() - 1, k + 1, V.getColumnCount() - 1, temp);
            }
        }

/*
      Scale the bidiagonal matrix so that its elements are
      real.
*/

        for (k = 0; k <= m; k++) {
            kk = k;
            xkk = X.unsafeGet(k, k);
            axkk = xkk.abs();
            X.unsafeSet(k, k, Complex.fromReal(axkk));
            d[kk] = axkk;
            scale = xkk.conj().divide(axkk);
            if (k < X.getColumnCount() - 1) {
                xkk1 = X.get(k, k + 1);
                xkk1.multiplyEquals(scale);
                X.unsafeSet(k, k + 1, xkk1);
            }
            scale = scale.conj();
            for (i = 0; i <= U.getRowCount() - 1; i++) {
                zr = U.unsafeGet(i, k).multiply(scale);
                U.unsafeSet(i, k, zr);
            }

            if (k < X.getColumnCount() - 1) {

                xkk1 = X.get(k, k + 1);
                axkk1 = xkk1.abs();
                X.unsafeSet(k, k + 1, Complex.fromReal(axkk1));
                e[kk] = axkk1;
                scale = xkk1.conj().divide(axkk1);
                if (k < X.getRowCount() - 1) {
                    xk1k1 = X.get(k + 1, k + 1);
                    xk1k1 = scale.multiply(xk1k1);
                    X.unsafeSet(k + 1, k + 1, xk1k1);
                }
                for (i = 0; i <= V.getRowCount() - 1; i++) {
                    zr = V.unsafeGet(i, k + 1).multiply(scale);
                    V.unsafeSet(i, k + 1, zr);
                }
            }
        }
        
/*
      If X has more columns than rows, rotate out the extra
      superdiagonal element.
*/
        if (X.getRowCount() < X.getColumnCount()) {
            t = e[m];
            for (k = m; k >= 0; k--) {
                ComplexPlaneRotationDense.genr(d[k], t, P);
                d[k] = P.zr;
                if (k != 0) {
                    t = P.sr * e[k - 1];
                    e[k - 1] = P.c * e[k - 1];
                }
                ComplexPlaneRotationDense.ap(V, P, 0, V.getRowCount() - 1, k, X.getRowCount() - 1 + 1);
                ComplexPlaneRotationDense.ap(X, P, 0, X.getRowCount() - 1, k, X.getRowCount() - 1 + 1);
            }
        }
/*
      Caculate the singular values of the bidiagonal matrix.
*/
        iu = m;
        iter = 0;
        while (true) {
/*
         These two loops determine the rows (il to iu) to
         iterate on.
*/
            while (iu > 0) {
                if (Math.abs(e[iu - 1]) >
                        1.0e-16 * (Math.abs(d[iu]) + Math.abs(d[iu - 1])))
                    break;
                e[iu - 1] = 0.;
                iter = 0;
                iu = iu - 1;
            }
            iter = iter + 1;
            if (iter > MAXITER) {
                throw new MaximumNumberOfIterationsExceededException("Maximum number of iterations exceeded.");
            }
            if (iu == 0) break;

            il = iu - 1;
            while (il > 0) {
                if (Math.abs(e[il - 1]) <=
                        1.0e-16 * (Math.abs(d[il]) + Math.abs(d[il - 1])))
                    break;
                il = il - 1;
            }
            if (il != 0) {
                e[il - 1] = 0.;
            }
/*
         Compute the shift (formulas adapted from LAPACK).
*/
            dmax = Math.max(Math.abs(d[iu]), Math.abs(d[iu - 1]));
            dmin = Math.min(Math.abs(d[iu]), Math.abs(d[iu - 1]));
            ea = Math.abs(e[iu - 1]);
            if (dmin == 0.) {
                shift = 0.;
            } else if (ea < dmax) {
                as = 1. + dmin / dmax;
                at = (dmax - dmin) / dmax;
                au = ea / dmax;
                au = au * au;
                shift = dmin * (2. / (Math.sqrt(as * as + au) + Math.sqrt(at * at + au)));
            } else {
                au = dmax / ea;
                if (au == 0.) {
                    shift = (dmin * dmax) / ea;
                } else {
                    as = 1. + dmin / dmax;
                    at = (dmax - dmin) / dmax;
                    t = 1. / (Math.sqrt(1. + (as * au) * (as * au)) +
                            Math.sqrt(1. + (at * au) * (at * au)));
                    shift = (t * dmin) * au;
                }
            }
/*
        Perform the implicitly shifted QR step.
*/
            t = Math.max(Math.max(Math.abs(d[il]), Math.abs(e[il])), shift);
            ds = d[il] / t;
            es = e[il] / t;
            ss = shift / t;
            ComplexPlaneRotationDense.genr((ds - ss) * (ds + ss), ds * es, P);
            for (i = il; i < iu; i++) {
                t = P.c * d[i] - P.sr * e[i];
                e[i] = P.sr * d[i] + P.c * e[i];
                d[i] = t;
                t = -P.sr * d[i + 1];
                d[i + 1] = P.c * d[i + 1];
                ComplexPlaneRotationDense.ap(V, P, 0, V.getRowCount() - 1, i, i + 1);
                ComplexPlaneRotationDense.genc(d[i], t, P);
                d[i] = P.zr;
                t = P.c * e[i] + P.sr * d[i + 1];
                d[i + 1] = P.c * d[i + 1] - P.sr * e[i];
                e[i] = t;
                ComplexPlaneRotationDense.aph(U, P, 0, U.getRowCount() - 1, i, i + 1);
                if (i != iu - 1) {
                    t = P.sr * e[i + 1];
                    e[i + 1] = P.c * e[i + 1];
                    ComplexPlaneRotationDense.genr(e[i], t, P);
                    e[i] = P.zr;
                }
            }
        }

/*
      Sort the singular values, setting negative values of d
      to positive.
*/
        for (k = m; k >= 0; k--) {
            if (d[k] < 0) {
                d[k] = -d[k];
                for (i = 0; i < X.getColumnCount(); i++) {
                    V.unsafeSet(i, k, V.unsafeGet(i, k).uminus());
                }
            }
            for (j = k; j < m; j++) {
                if (d[j] < d[j + 1]) {
                    t = d[j];
                    d[j] = d[j + 1];
                    d[j + 1] = t;
                    for (i = 0; i < X.getRowCount(); i++) {
                        t = U.unsafeGet(i, j).real();
                        double a = U.unsafeGet(i, j + 1).real();
                        double b = t;
                        t = U.unsafeGet(i, j).imag();
                        double c = U.unsafeGet(i, j + 1).imag();
                        U.unsafeSet(i, j, new Complex(a, c));
                        U.unsafeSet(i, j + 1, new Complex(b, t));
                    }
                    for (i = 0; i < X.getColumnCount(); i++) {
                        t = V.unsafeGet(i, j).real();
                        double a = V.unsafeGet(i, j + 1).real();
                        double b = t;
                        t = V.unsafeGet(i, j).imag();
                        double c = V.unsafeGet(i, j + 1).imag();
                        V.unsafeSet(i, j, new Complex(a, c));
                        V.unsafeSet(i, j + 1, new Complex(b, t));
                    }
                }
            }
        }
/*
      Return the decompostion;
*/
        for(i = 0; i < mc; i++) {
            S.unsafeSet(i, i, d[i]);
        }
    }

    public ComplexMatrixDense getU() {
        return U;
    }

    public MatrixDense getS() {
        return S;
    }

    public double[] getSingularValues() {
        return d;
    }

    public ComplexMatrixDense getV() {
        return V;
    }

    /**
     * Effective numerical matrix rank
     *
     * @return Number of nonnegligible singular values.
     */

    public int rank() {
        double eps = ConstantsETK.DOUBLE_EPS;
        double[] values = S.diag();
        double tol = Math.max(rows, cols) * values[0] * eps;
        int r = 0;
        for (double v : values) {
            if (v > tol) {
                r++;
            }
        }
        return r;
    }
}
