package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;

public class ComplexHouseholderTransformations {

    /**
     * Generates a Householder transformation from within the part of
     * column c of a ComplexMatrixDense (altered) extending from rows
     * r1 to r2.  The method overwrites the
     * column with the result of applying the transformation.
     *
     * @param A  The matrix from which the transformation
     *           is to be generated (altered)
     * @param r1 The index of the row in which the generating column
     *           begins
     * @param r2 The index of the row in which the generating column
     *           ends
     * @param c  The index of the generating column
     * @return A Complex[] of length r2-r1+1
     * containing the Householder vector
     * @throws RuntimeException Passed from below.
     */
    public static Complex[] genc(ComplexMatrixDense A, int r1, int r2, int c)
            throws RuntimeException {

        int i, ru;
        double norm;
        double s;
        Complex scale;
        Complex t;
        Complex t1;

        ru = r2 - r1 + 1;

        Complex[] u = new Complex[r2 - r1 + 1];

        for (i = r1; i <= r2; i++) {
            u[i - r1] = new Complex(A.unsafeGet(i, c));
            A.unsafeSet(i, c, new Complex());
        }

        norm = ComplexArrays.normFro(u);

        if (r1 == r2 || norm == 0) {
            A.unsafeSet(r1, c, u[0].uminus());
            u[0] =  Complex.fromReal(Math.sqrt(2));
            return u;
        }

        scale = new Complex(1 / norm, 0);

        if (u[0].real() != 0 || u[0].imag() != 0) {
            t = u[0].copy();
            t1 = t.conj();
            t = t1.divide(t.abs());
            scale.multiplyEquals(t);
        }

        t = Complex.fromReal(1).divide(scale).uminus();;
        A.unsafeSet(r1, c, t);

        for (i = 0; i < ru; i++) {
            u[i].multiplyEquals(scale);
        }

        u[0] = new Complex(u[0].real() + 1, 0);
        s = Math.sqrt(1 / u[0].real());

        for (i = 0; i < ru; i++) {
            u[i] = new Complex(s * u[i].real(), s * u[i].imag());
        }

        return u;
    }

    /**
     * Generates a Householder transformation from within the part of row
     * r of a ComplexMatrixDense (altered) extending from columns c1 to
     * c2.  The method overwrites the row with the result
     * of applying the transformation.
     *
     * @param A  The matrix from which the transformation
     *           is to be generated (altered)
     * @param r  The index of the generating row
     * @param c1 The index of the column in which the generating row
     *           begins
     * @param c2 The index of the column in which the generating row
     *           ends
     * @return A Complex[] of length r2-r1+1
     * containing the Householder vector
     * @throws RuntimeException Passed from below.
     */
    public static Complex[] genr(ComplexMatrixDense A, int r, int c1, int c2)
            throws RuntimeException {

        int j, cu;
        double norm, s;
        Complex scale;
        Complex t = new Complex();
        Complex t1 = new Complex();

        cu = c2 - c1 + 1;

        Complex[] u = new Complex[cu];

        for (j = c1; j <= c2; j++) {
            u[j - c1] = new Complex(A.unsafeGet(r, j));
            A.unsafeSet(r, j, new Complex());
        }

        norm = ComplexArrays.normFro(u);

        if (c1 == c2 || norm == 0) {
            A.unsafeSet(r, c1, u[0].uminus());
            u[0] = Complex.fromReal(Math.sqrt(2));
            return u;
        }

        scale = new Complex(1 / norm, 0);

        if (u[0].real() != 0 || u[0].imag() != 0) {
            t = u[0];
            // TODO
            //scale.multiplyEquals();
            //scale.Times(scale, t.Div(t1.Conj(t), Complex.abs(t)));
        }


        //A.unsafeSet(r, c1, t.Minus(t.Div(Complex.ONE, scale)));

        for (j = 0; j < cu; j++) {
            //u.Times(j, scale);
        }

        u[0] = Complex.fromReal(u[0].real() + 1);
        s = Math.sqrt(1 / u[0].real());

        for (j = 0; j < cu; j++) {
            u[j] = new Complex(s * u[j].real(), -s * u[j].imag());
        }

        return u;
    }

    /**
     * Premultiplies the Householder transformation contained in a
     * Complex[] into a ComplexMatrixDense A[r1:r2,c1:c2] and overwrites
     * ComplexMatrixDense A[r1:r2,c1:c2] with the results.  If r1 &gt; r2
     * or c1 &gt; c2 the method does nothing.
     *
     * @param u  The Householder vector
     * @param A  The ComplexMatrixDense to which the transformation
     *           is to be applied (altered)
     * @param r1 The index of the first row to which the transformation
     *           is to be applied
     * @param r2 The index of the last row to which the transformation
     *           is to be applied
     * @param c1 The index of the first column to which the transformation
     *           is index of the to be applied
     * @param c2 The index of the last column to which the transformation
     *           is to be applied
     * @param v  A work array of length at least c2-c1+1
     * @return The transformed ComplexMatrixDense A
     * @throws RuntimeException Thrown if either u or v is too short.
     */
    public static ComplexMatrixDense ua(Complex[] u, ComplexMatrixDense A, int r1, int r2, int c1, int c2, Complex[] v)
            throws RuntimeException {

        int i, j;


        if (r2 < r1 || c2 < c1) {
            return A;
        }

        if (r2 - r1 + 1 > u.length) {
            throw new RuntimeException
                    ("Householder vector too short.");
        }

        if (c2 - c1 + 1 > v.length) {
            throw new RuntimeException
                    ("Work vector too short.");
        }

        for (j = c1; j <= c2; j++) {
            v[j - c1] = new Complex();
        }

        for (i = r1; i <= r2; i++) {
            for (j = c1; j <= c2; j++) {
                double a = v[j - c1].real();
                double b = u[i - r1].real();
                double c = A.unsafeGet(i, j).real();
                double d = u[i - r1].imag();
                double e = A.unsafeGet(i, j).imag();

                double aa = v[j - c1].imag();
                double bb = u[i - r1].real();
                double cc = A.unsafeGet(i, j).imag();
                double dd = u[i - r1].imag();
                double ee = A.unsafeGet(i, j).real();
                v[j - c1] = new Complex(a + b * c + d * e, aa + bb * cc - dd * ee);
            }
        }

        for (i = r1; i <= r2; i++) {
            for (j = c1; j <= c2; j++) {
                double a = A.unsafeGet(i, j).real();
                double b = u[i - r1].real();
                double c = v[j - c1].real();
                double d = u[i - r1].imag();
                double e = v[j - c1].imag();

                double aa = A.unsafeGet(i, j).imag();
                double bb = u[i - r1].real();
                double cc = v[j - c1].imag();
                double dd = u[i - r1].imag();
                double ee= v[j - c1].real();
                A.unsafeSet(i, j, new Complex(a - b * c + d * e, aa - bb * cc - dd * ee));
            }
        }
        return A;
    }

    /**
     * Premultiplies the Householder transformation contained in a
     * Complex[] into a ComplexMatrixDense A[r1:r2,c1:c2] and overwrites
     * ComplexMatrixDense A[r1:r2,c1:c2] with the results.  If r1 &gt; r2
     * or c1 &gt; c2 the method does nothing.
     *
     * @param u  The Householder vector
     * @param A  The ComplexMatrixDense to which the transformation
     *           is to be applied (altered)
     * @param r1 The index of the first row to which the transformation
     *           is to be applied
     * @param r2 The index of the last row to which the transformation
     *           is to be applied
     * @param c1 The index of the first column to which the transformation
     *           is index of the to be applied
     * @param c2 The index of the last column to which the transformation
     *           is to be applied
     * @return The transformed ComplexMatrixDense A
     * @throws RuntimeException Passed from below.
     */
    public static ComplexMatrixDense ua(Complex[] u, ComplexMatrixDense A, int r1, int r2, int c1, int c2)
            throws RuntimeException {

        if (c1 > c2) {
            return A;
        }

        return ua(u, A, r1, r2, c1, c2, new Complex[c2 - c1 + 1]);
    }


    /**
     * Postmultiplies the Householder transformation contained in a
     * Complex[] into a ComplexMatrixDense A[r1:r2,c1:c2] and overwrites
     * ComplexMatrixDense A[r1:r2,c1:c2] with the results.  If r1 &gt; r2
     * or c1 &gt; c2 the method does nothing.
     *
     * @param u  The Householder vector
     * @param A  The ComplexMatrixDense to which the transformation
     *           is to be applied (altered)
     * @param r1 The index of the first row to which the transformation
     *           is to be applied
     * @param r2 The index of the last row to which the transformation
     *           is to be applied
     * @param c1 The index of the first column to which the transformation
     *           is index of the to be applied
     * @param c2 The index of the last column to which the transformation
     *           is to be applied
     * @param v  A work array of length at least c2-c1+1
     * @return The transformed ComplexMatrixDense A
     * @throws RuntimeException Thrown if either u or v is too short.
     */
    public static ComplexMatrixDense au(ComplexMatrixDense A, Complex[] u, int r1, int r2, int c1, int c2, Complex[] v)
            throws RuntimeException {

        int i, j, cu;

        if (r2 < r1 || c2 < c1) {
            return A;
        }

        if (c2 - c1 + 1 > u.length) {
            throw new RuntimeException
                    ("Householder vector too short.");
        }

        if (r2 - r1 + 1 > v.length) {
            throw new RuntimeException
                    ("Work vector too short.");
        }

        for (i = r1; i <= r2; i++) {
            v[i - r1] = new Complex();
            for (j = c1; j <= c2; j++) {
                double a = v[i - r1].real();
                double b = A.unsafeGet(i, j).real();
                double c = u[j - c1].real();
                double d = A.unsafeGet(i, j).imag();
                double e = u[j - c1].imag();

                double aa = v[i - r1].imag();
                double bb = A.unsafeGet(i, j).real();
                double cc = u[j - c1].imag();
                double dd = A.unsafeGet(i, j).imag();
                double ee = u[j - c1].real();
                v[i - r1] = new Complex(a + b * c - d * e, aa + bb * cc + dd * ee);
            }
        }
        for (i = r1; i <= r2; i++) {
            for (j = c1; j <= c2; j++) {
                double a = A.unsafeGet(i, j).real();
                double b = v[i - r1].real();
                double c = u[j - c1].real();
                double d = v[i - r1].imag();
                double e = u[j - c1].imag();

                double aa = A.unsafeGet(i, j).imag();
                double bb = v[i - r1].real();
                double cc = u[j - c1].imag();
                double dd = v[i - r1].imag();
                double ee = u[j - c1].real();
                A.unsafeSet(i, j, new Complex(a - b * c - d * e, aa + bb * cc - dd * ee));
            }
        }
        return A;
    }


    /**
     * Postmultiplies the Householder transformation contained in a
     * Complex[] into a ComplexMatrixDense A[r1:r2,c1:c2] and overwrites
     * ComplexMatrixDense A[r1:r2,c1:c2] with the results.  If r1 &gt; r2
     * or c1 &gt; c2 the method does nothing.
     *
     * @param u  The Householder vector
     * @param A  The ComplexMatrixDense to which the transformation
     *           is to be applied (altered)
     * @param r1 The index of the first row to which the transformation
     *           is to be applied
     * @param r2 The index of the last row to which the transformation
     *           is to be applied
     * @param c1 The index of the first column to which the transformation
     *           is index of the to be applied
     * @param c2 The index of the last column to which the transformation
     *           is to be applied
     * @return The transformed ComplexMatrixDense A
     * @throws RuntimeException Passed from below.
     */

    public static ComplexMatrixDense au(ComplexMatrixDense A, Complex[] u, int r1, int r2, int c1, int c2)
            throws RuntimeException {

        if (r2 < r1) {
            return A;
        }

        return au(A, u, r1, r2, c1, c2, new Complex[r2 - r1 + 1]);
    }

}
