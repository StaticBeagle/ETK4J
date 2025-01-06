package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

/**
 * Rot generates and manipulates plane rotations.  Given a 2-vector
 * compontents are x and y, there is a unitary matrix P such that
 * <pre>
 *      P|x| =  |   c      s||x| = |z|
 *       |y|    |-conj(s)  c||y|   |0|
 * </pre>
 * The number c, which is always real, is the cosine of the rotation.
 * The number s, which may be complex is the sine of the rotation.
 * <p>
 * Comments: This suite will eventually contain methods for real
 * rotations (two are already in place).  The only difference
 * between real and complex rotations is that si and zi are zero
 * for the former.  The final routines will do the efficient thing.
 *
 * @author G. W. Stewart
 * @version Pre-alpha
 */
public class ComplexPlaneRotationDense {

    /**
     * The cosine of the rotation
     */
    protected double c;
    /**
     * The real part of the sine of the rotation
     */
    public double sr;
    /**
     * The imaginary part of the sine of the rotation
     */
    public double si;
    /**
     * The real part of the first component of the transformed vector
     */
    public double zr;
    /**
     * The imaginary part of the first component of the transformed vector
     */
    public double zi;

    /* // TODO fix the javadocs
     * Given the real and imaginary parts of a 2-vector, genc returns
     * a plane rotation P such that
     * <pre>
     *      P|x| =  |   c      s||x| = |z|
     *       |y|    |-conj(s)  c||y|   |0|
     * </pre>
     *
     * @param xr The real part of the first component of the 2-vector
     * @param xi The imaginary part of the first component of the 2-vector
     * @param yr The real part of the second component of the 2-vector
     * @param yi The imaginary part of the second component of the 2-vector
     * @return The rotation
     */

    public static ComplexPlaneRotationDense genc(Complex x, Complex y) {

        double s, absx, absxy;

        ComplexPlaneRotationDense P = new ComplexPlaneRotationDense();

        if (x.real() == 0 && x.imag() == 0) {
            P.c = 0.;
            P.sr = 1.;
            P.si = 0.;
            P.zr = y.real();
            P.zi = y.imag();
            return P;
        }
        double xr = x.real();
        double xi = x.imag();

        double yr = y.real();
        double yi = y.imag();

        s = Math.abs(xr) + Math.abs(xi);
        absx = s * Math.sqrt((xr / s) * (xr / s) + (xi / s) * (xi / s));
        s = Math.abs(s) + Math.abs(yr) + Math.abs(yi);
        absxy = s * Math.sqrt((absx / s) * (absx / s) + (yr / s) * (yr / s) + (yi / s) * (yi / s));
        P.c = absx / absxy;
        xr = xr / absx;
        xi = xi / absx;
        P.sr = (xr * yr + xi * yi) / absxy;
        P.si = (xi * yr - xr * yi) / absxy;
        P.zr = xr * absxy;
        P.zi = xi * absxy;
        return P;
    }

    /* TODO
     * Given the real and imaginary parts of a 2-vector, genc generates
     * a plane rotation P such that
     * <pre>
     *      P|x| =  |   c      s||x| = |z|
     *       |y|    |-conj(s)  c||y|   |0|
     * </pre>
     *
     * @param xr The real part of the first component of the 2-vector
     * @param xi The imaginary part of the first component of the 2-vector
     * @param yr The real part of the second component of the 2-vector
     * @param yi The imaginary part of the second component of the 2-vector
     * @param P  The rotation (must be initialized)
     */

    public static void genc(Complex x, Complex y, ComplexPlaneRotationDense P) {

        double s, absx, absxy;

        if (x.real() == 0 && x.imag() == 0) {
            P.c = 0.;
            P.sr = 1.;
            P.si = 0.;
            P.zr = y.real();
            P.zi = y.imag();
            return;
        }
        double xr = x.real();
        double xi = x.imag();
        double yr = y.real();
        double yi = y.imag();
        s = Math.abs(xr) + Math.abs(xi);
        absx = s * Math.sqrt((xr / s) * (xr / s) + (xi / s) * (xi / s));
        s = Math.abs(s) + Math.abs(yr) + Math.abs(yi);
        absxy = s * Math.sqrt((absx / s) * (absx / s) + (yr / s) * (yr / s) + (yi / s) * (yi / s));
        P.c = absx / absxy;
        xr = xr / absx;
        xi = xi / absx;
        P.sr = (xr * yr + xi * yi) / absxy;
        P.si = (xi * yr - xr * yi) / absxy;
        P.zr = xr * absxy;
        P.zi = xi * absxy;
    }

    /**
     * Given a real 2-vector, genc returns
     * a real plane rotation P such that
     * <pre>
     *      P|x| =  | c  s||x| = |z|
     *       |y|    |-s  c||y|   |0|
     * </pre>
     *
     * @param x The first component of the two vector
     * @param y The second component of the two vector
     * @return The rotation
     */
    public static ComplexPlaneRotationDense genc(double x, double y) {

        ComplexPlaneRotationDense P = new ComplexPlaneRotationDense();

        P.si = 0.;
        P.zi = 0.;

        if (x == 0 & y == 0) {
            P.c = 1;
            P.sr = 0.;
            P.zr = 0.;
            return P;
        }

        double s = Math.abs(x) + Math.abs(y);
        P.zr = s * Math.sqrt((x / s) * (x / s) + (y / s) * (y / s));
        P.c = x / P.zr;
        P.sr = y / P.zr;
        return P;
    }

    /**
     * Given a real 2-vectc, genc generates
     * a real plane rotation P such that
     * <pre>
     *      P|x| =  | c  s||x| = |z|
     *       |y|    |-s  c||y|   |0|
     * </pre>
     *
     * @param x The first component of the two vector
     * @param y The second component of the two vector
     * @param P The plane rotation
     */
    public static void genc(double x, double y, ComplexPlaneRotationDense P) {

        P.si = 0.;
        P.zi = 0.;

        if (x == 0 & y == 0) {
            P.c = 1;
            P.sr = 0.;
            P.zr = 0.;
            return;
        }

        double s = Math.abs(x) + Math.abs(y);
        P.zr = s * Math.sqrt((x / s) * (x / s) + (y / s) * (y / s));
        P.c = x / P.zr;
        P.sr = y / P.zr;
    }

    /**
     * Given a Zmat A, genc returns a plane rotation that on
     * premultiplication into rows ii1 and ii2
     * annihilates A(ii2,jj).  The element A(ii2,jj) is
     * overwriten by zero and the element A(ii1,jj) is
     * overwritten by its transformed value.
     *
     * @param A   The Zmat (altered)
     * @param ii1 The row index of the first element
     * @param ii2 The row index of the second element (the
     *            one that is annihilated
     * @param jj  The column index of the elements
     * @return The plane rotation
     */
    public static ComplexPlaneRotationDense genc(ComplexMatrixDense A, int ii1, int ii2, int jj) {

        int i1 = ii1;
        int i2 = ii2;
        int j = jj;

        ComplexPlaneRotationDense P = ComplexPlaneRotationDense.genc(A.unsafeGet(i1, j), A.unsafeGet(i2, j));
        A.unsafeSet(i1, j, new Complex(P.zr, P.zi));
        A.unsafeSet(i2, j, new Complex());
        return P;
    }

    /**
     * Given a Zmat A, genc generates a plane rotation that on
     * premultiplication into rows ii1 and ii2
     * annihilates A(ii2,jj).  The element A(ii2,jj) is
     * overwriten by zero and the element A(ii1,jj) is
     * overwritten by its transformed value.
     *
     * @param A   The Zmat (altered)
     * @param ii1 The row index of the first element
     * @param ii2 The row index of the second element (the
     *            one that is annihilated
     * @param jj  The column index of the elements
     * @param P   The plane rotation (must be initialized)
     */
    public static void genc(ComplexMatrixDense A, int ii1, int ii2, int jj, ComplexPlaneRotationDense P) {

        int i1 = ii1;
        int i2 = ii2;
        int j = jj;

        ComplexPlaneRotationDense.genc(A.unsafeGet(i1, j), A.unsafeGet(i2, j), P);
        A.unsafeSet(i1, j, new Complex(P.zr, P.zi));
        A.unsafeSet(i2, j, new Complex());
    }

    /* TODO
     * Given the real and imaginary parts of a 2-vector, genr returns
     * a plane rotation such that
     * <pre>
     *      |x y|P = |x y||   c      s||x| = |z 0|
     *                    |-conj(s)  c||y|
     * </pre>
     *
     * @param xr The real part of the first component of the 2-vector
     * @param xi The imaginary part of the first component of the 2-vector
     * @param yr The real part of the second component of the 2-vector
     * @param yi The imaginary part of the second component of the 2-vector
     * @return The rotation
     */

    public static ComplexPlaneRotationDense genr(Complex x, Complex y) {

        double s, absx, absxy;

        ComplexPlaneRotationDense P = new ComplexPlaneRotationDense();

        if (x.real() == 0 && x.imag() == 0) {
            P.c = 0.;
            P.sr = 1.;
            P.si = 0.;
            P.zr = y.real();
            P.zi = y.imag();
            return P;
        }
        double xr = x.real();
        double xi = x.imag();
        double yr = y.real();
        double yi = y.imag();

        s = Math.abs(xr) + Math.abs(xi);
        absx = s * Math.sqrt((xr / s) * (xr / s) + (xi / s) * (xi / s));
        s = Math.abs(s) + Math.abs(yr) + Math.abs(yi);
        absxy = s * Math.sqrt((absx / s) * (absx / s) + (yr / s) * (yr / s) + (yi / s) * (yi / s));
        P.c = absx / absxy;
        xr = xr / absx;
        xi = xi / absx;
        P.sr = -(xr * yr + xi * yi) / absxy;
        P.si = (xi * yr - xr * yi) / absxy;
        P.zr = xr * absxy;
        P.zi = xi * absxy;
        return P;
    }

    /* TODO
     * Given the real and imaginary parts of a 2-vector, genr generates
     * a plane rotation such that
     * <pre>
     *      |x y|P = |x y||   c      s||x| = |z 0|
     *                    |-conj(s)  c||y|
     * </pre>
     *
     * @param xr The real part of the first component of the 2-vector
     * @param xi The imaginary part of the first component of the 2-vector
     * @param yr The real part of the second component of the 2-vector
     * @param yi The imaginary part of the second component of the 2-vector
     * @param P  The plane rotation (must be initialized)
     */

    public static void genr(Complex x, Complex y, ComplexPlaneRotationDense P) {

        double s, absx, absxy;

        if (x.real() == 0 && x.imag() == 0) {
            P.c = 0.;
            P.sr = 1.;
            P.si = 0.;
            P.zr = y.real();
            P.zi = y.imag();
            return;
        }
        double xr = x.real();
        double xi = x.imag();
        double yr = y.real();
        double yi = y.imag();

        s = Math.abs(xr) + Math.abs(xi);
        absx = s * Math.sqrt((xr / s) * (xr / s) + (xi / s) * (xi / s));
        s = Math.abs(s) + Math.abs(yr) + Math.abs(yi);
        absxy = s * Math.sqrt((absx / s) * (absx / s) + (yr / s) * (yr / s) + (yi / s) * (yi / s));
        P.c = absx / absxy;
        xr = xr / absx;
        xi = xi / absx;
        P.sr = -(xr * yr + xi * yi) / absxy;
        P.si = (xi * yr - xr * yi) / absxy;
        P.zr = xr * absxy;
        P.zi = xi * absxy;
    }

    /**
     * Given a Zmat A, genr returns a plane rotation that on
     * postmultiplication into column jj1 and jj2
     * annihilates A(ii,jj2).  The element A(ii,jj2) is
     * overwirten by zero and the element A(ii,jj1) is
     * overwritten by its transformed value.
     *
     * @param A   The Zmat (altered)
     * @param ii  The index of the row containing the elements
     * @param jj1 The column index of the first element
     * @param jj2 The column index of the second element (the
     *            one that is annihilated)
     * @return The rotation
     */

    public static ComplexPlaneRotationDense genr(ComplexMatrixDense A, int ii, int jj1, int jj2) {

        int i = ii;
        int j1 = jj1;
        int j2 = jj2;

        ComplexPlaneRotationDense P = ComplexPlaneRotationDense.genr(A.unsafeGet(i, j1), A.unsafeGet(i, j2));
        A.unsafeSet(i, j1, new Complex(P.zr, P.zi));
        A.unsafeSet(i, j2, new Complex());
        return P;
    }

    /**
     * Given a Zmat A, genr generates a plane rotation that on
     * postmultiplication into column jj1 and jj2
     * annihilates A(ii,jj2).  The element A(ii,jj2) is
     * overwirten by zero and the element A(ii,jj1) is
     * overwritten by its transformed value.
     *
     * @param A   The Zmat (altered)
     * @param ii  The index of the row containing the elements
     * @param jj1 The column index of the first element
     * @param jj2 The column index of the second element (the
     *            one that is annihilated)
     * @param P   The rotation
     */

    public static void genr(ComplexMatrixDense A, int ii, int jj1, int jj2, ComplexPlaneRotationDense P) {
        ComplexPlaneRotationDense.genr(A.unsafeGet(ii, jj1), A.unsafeGet(ii, jj2), P);
        A.unsafeSet(ii, jj1, new Complex(P.zr, P.zi));
        A.unsafeSet(ii, jj2, new Complex());
        A.unsafeSet(ii, jj2, new Complex());
    }

    /**
     * Given  a real 2-vector, genr returns
     * a plane rotation such that
     * <pre>
     *      |x y|P = |x y|| c  s||x| = |z 0|
     *                    |-s  c||y|
     * </pre>
     *
     * @param x The first component of the 2-vector
     * @param y The second component of the 2-vector
     * @return The rotation
     */
    public static ComplexPlaneRotationDense genr(double x, double y) {

        ComplexPlaneRotationDense P = new ComplexPlaneRotationDense();

        P.si = 0.;
        P.zi = 0.;

        double s = Math.abs(x) + Math.abs(y);

        if (s == 0.) {
            P.c = 1.;
            P.sr = 0.;
            P.zr = 0.;
            return P;
        }

        P.zr = s * Math.sqrt((x / s) * (x / s) + (y / s) * (y / s));
        P.c = x / P.zr;
        P.sr = -y / P.zr;
        return P;
    }


    /**
     * Given  a real 2-vector, genr generates
     * a plane rotation such that
     * <pre>
     *      |x y|P = |x y|| c  s||x| = |z 0|
     *                    |-s  c||y|
     * </pre>
     *
     * @param x The first component of the 2-vector
     * @param y The second component of the 2-vector
     * @param P The rotation
     */
    public static void genr(double x, double y, ComplexPlaneRotationDense P) {

        P.si = 0.;
        P.zi = 0.;

        double s = Math.abs(x) + Math.abs(y);

        if (s == 0.) {
            P.c = 1.;
            P.sr = 0.;
            P.zr = 0.;
            return;
        }

        P.zr = s * Math.sqrt((x / s) * (x / s) + (y / s) * (y / s));
        P.c = x / P.zr;
        P.sr = -y / P.zr;
    }


    /**
     * Multiplies rows (ii1,jj1:jj2) and (ii2,jj1:jj2)
     * of a Zmat (altered) by a plane rotation.
     *
     * @param P   The plane rotation
     * @param A   The Zmat (altered)
     * @param ii1 The row index of the first row.
     * @param ii2 The row index of the second row.
     * @param jj1 The first index of the range of the rows
     * @param jj2 The second index of the range of the rows
     */


    public static void pa(ComplexPlaneRotationDense P, ComplexMatrixDense A, int ii1, int ii2, int jj1, int jj2) {

        double t1r, t1i, t2r, t2i;

        for (int j = jj1; j <= jj2; j++) {
            t1r = P.c * A.unsafeGet(ii1, j).real() + P.sr * A.unsafeGet(ii2, j).real() - P.si * A.unsafeGet(ii2, j).imag();
            t1i = P.c * A.unsafeGet(ii1, j).imag() + P.sr * A.unsafeGet(ii2, j).imag() + P.si * A.unsafeGet(ii2, j).real();
            t2r = P.c * A.unsafeGet(ii2, j).real() - P.sr * A.unsafeGet(ii1, j).real() - P.si * A.unsafeGet(ii1, j).imag();
            t2i = P.c * A.unsafeGet(ii2, j).imag() - P.sr * A.unsafeGet(ii1, j).imag() + P.si * A.unsafeGet(ii1, j).real();
            A.unsafeSet(ii1, j, new Complex(t1r, t1i));
            A.unsafeSet(ii2, j, new Complex(t2r, t2i));
        }
    }

    /**
     * Multiplies rows (ii1,jj1:jj2) and (ii2,jj1:jj2)
     * of a Zmat (altered) by the conjugate transpose of a plane rotation.
     *
     * @param P   The plane rotation
     * @param A   The Zmat (altered)
     * @param ii1 The row index of the first row.
     * @param ii2 The row index of the second row.
     * @param jj1 The first index of the range of the rows
     * @param jj2 The second index of the range of the rows
     */
    public static void pha(ComplexPlaneRotationDense P, ComplexMatrixDense A, int ii1, int ii2, int jj1, int jj2) {

        double t1r, t1i, t2r, t2i;

        int i1 = ii1;
        int i2 = ii2;
        int j1 = jj1;
        int j2 = jj2;


        for (int j = j1; j <= j2; j++) {
            t1r = P.c * A.unsafeGet(i1, j).real() - P.sr * A.unsafeGet(i2, j).real() + P.si * A.unsafeGet(i2, j).imag();
            t1i = P.c * A.unsafeGet(i1, j).imag() - P.sr * A.unsafeGet(i2, j).imag() - P.si * A.unsafeGet(i2, j).real();
            t2r = P.c * A.unsafeGet(i2, j).real() + P.sr * A.unsafeGet(i1, j).real() + P.si * A.unsafeGet(i1, j).imag();
            t2i = P.c * A.unsafeGet(i2, j).imag() + P.sr * A.unsafeGet(i1, j).imag() - P.si * A.unsafeGet(i1, j).real();
            A.unsafeSet(i1, j, new Complex(t1r, t1i));
            A.unsafeSet(i2, j, new Complex(t2r, t2i));
        }
    }


    /**
     * Multiplies columns (ii1:ii2,jj1) and A(ii2:ii2,jj1)
     * of a Zmat (altered) by a plane rotation.
     *
     * @param A   The Zmat (altered)
     * @param P   The rotation
     * @param ii1 The first index of the column range
     * @param ii2 The second index of the column range
     * @param jj1 The index of the first column
     * @param jj2 The index of the second column
     */

    public static void ap(ComplexMatrixDense A, ComplexPlaneRotationDense P, int ii1, int ii2, int jj1, int jj2) {

        double t1r, t1i, t2r, t2i;
        for (int i = ii1; i <= ii2; i++) {
            t1r = P.c * A.unsafeGet(i, jj1).real() - P.sr * A.unsafeGet(i, jj2).real() - P.si * A.unsafeGet(i, jj2).imag();
            t1i = P.c * A.unsafeGet(i, jj1).imag() - P.sr * A.unsafeGet(i, jj2).imag() + P.si * A.unsafeGet(i, jj2).real();
            t2r = P.c * A.unsafeGet(i, jj2).real() + P.sr * A.unsafeGet(i, jj1).real() - P.si * A.unsafeGet(i, jj1).imag();
            t2i = P.c * A.unsafeGet(i, jj2).imag() + P.sr * A.unsafeGet(i, jj1).imag() + P.si * A.unsafeGet(i, jj1).real();
            A.unsafeSet(i, jj1, new Complex(t1r, t1i));
            A.unsafeSet(i, jj2, new Complex(t2r, t2i));
        }
    }

    /**
     * Multiplies columns (ii1:ii2,jj1) and A(ii2:ii2,jj1)
     * of a Zmat (altered) by the conjugate transpose of plane rotation.
     *
     * @param A   The Zmat (altered)
     * @param P   The rotation
     * @param ii1 The first index of the column range
     * @param ii2 The second index of the column range
     * @param jj1 The index of the first column
     * @param jj2 The index of the second column
     */

    public static void aph(ComplexMatrixDense A, ComplexPlaneRotationDense P, int ii1, int ii2, int jj1, int jj2) {

        double t1r, t1i, t2r, t2i;
        for (int i = ii1; i <= ii2; i++) {
            t1r = P.c * A.unsafeGet(i, jj1).real() + P.sr * A.unsafeGet(i, jj2).real() + P.si * A.unsafeGet(i, jj2).imag();
            t1i = P.c * A.unsafeGet(i, jj1).imag() + P.sr * A.unsafeGet(i, jj2).imag() - P.si * A.unsafeGet(i, jj2).real();
            t2r = P.c * A.unsafeGet(i, jj2).real() - P.sr * A.unsafeGet(i, jj1).real() + P.si * A.unsafeGet(i, jj1).imag();
            t2i = P.c * A.unsafeGet(i, jj2).imag() - P.sr * A.unsafeGet(i, jj1).imag() - P.si * A.unsafeGet(i, jj1).real();
            A.unsafeSet(i, jj1, new Complex(t1r, t1i));
            A.unsafeSet(i, jj2, new Complex(t2r, t2i));
        }
    }
}
