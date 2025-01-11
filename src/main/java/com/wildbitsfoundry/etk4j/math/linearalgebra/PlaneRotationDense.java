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
public class PlaneRotationDense {

    /**
     * The cosine of the rotation
     */
    protected double c;
    /**
     * The real part of the sine of the rotation
     */
    public double sr;
    /**
     * The real part of the first component of the transformed vector
     */
    public double zr;


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

    public static PlaneRotationDense genc(double x, double y) {

        double s, absx, absxy;

        PlaneRotationDense P = new PlaneRotationDense();

        if (x == 0) {
            P.c = 0.;
            P.sr = 1.;
            P.zr = y;
            return P;
        }
        double xr = x;

        s = Math.abs(xr);
        absx = s * Math.sqrt((xr / s) * (xr / s));
        s = Math.abs(s) + Math.abs(y);
        absxy = s * Math.sqrt((absx / s) * (absx / s) + (y / s) * (y / s));
        P.c = absx / absxy;
        xr = xr / absx;
        P.sr = (xr * y) / absxy;
        P.zr = xr * absxy;
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

    public static void genc(double x, double y, PlaneRotationDense P) {

        double s, absx, absxy;

        if (x == 0) {
            P.c = 0.;
            P.sr = 1.;
            P.zr = y;
            return;
        }
        double xr = x;
        double yr = y;
        s = Math.abs(xr);
        absx = s * Math.sqrt((xr / s) * (xr / s));
        s = Math.abs(s) + Math.abs(yr);
        absxy = s * Math.sqrt((absx / s) * (absx / s) + (yr / s) * (yr / s));
        P.c = absx / absxy;
        xr = xr / absx;
        P.sr = (xr * yr) / absxy;
        P.zr = xr * absxy;
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
//    public static PlaneRotationDense genc(double x, double y) {
//
//        PlaneRotationDense P = new PlaneRotationDense();
//
//        if (x == 0 & y == 0) {
//            P.c = 1;
//            P.sr = 0.;
//            P.zr = 0.;
//            return P;
//        }
//
//        double s = Math.abs(x) + Math.abs(y);
//        P.zr = s * Math.sqrt((x / s) * (x / s) + (y / s) * (y / s));
//        P.c = x / P.zr;
//        P.sr = y / P.zr;
//        return P;
//    }

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
//    public static void genc(double x, double y, PlaneRotationDense P) {
//
//        if (x == 0 & y == 0) {
//            P.c = 1;
//            P.sr = 0.;
//            P.zr = 0.;
//            return;
//        }
//
//        double s = Math.abs(x) + Math.abs(y);
//        P.zr = s * Math.sqrt((x / s) * (x / s) + (y / s) * (y / s));
//        P.c = x / P.zr;
//        P.sr = y / P.zr;
//    }

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
    public static PlaneRotationDense genc(MatrixDense A, int ii1, int ii2, int jj) {

        int i1 = ii1;
        int i2 = ii2;
        int j = jj;

        PlaneRotationDense P = PlaneRotationDense.genc(A.unsafeGet(i1, j), A.unsafeGet(i2, j));
        A.unsafeSet(i1, j, P.zr);
        A.unsafeSet(i2, j, 0);
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
    public static void genc(MatrixDense A, int ii1, int ii2, int jj, PlaneRotationDense P) {

        PlaneRotationDense.genc(A.unsafeGet(ii1, jj), A.unsafeGet(ii2, jj), P);
        A.unsafeSet(ii1, jj, P.zr);
        A.unsafeSet(ii2, jj, 0);
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

    public static PlaneRotationDense genr(double x, double y) {

        double s, absx, absxy;

        PlaneRotationDense P = new PlaneRotationDense();

        if (x == 0) {
            P.c = 0.;
            P.sr = 1.;
            P.zr = y;
            return P;
        }
        double xr = x;
        double yr = y;

        s = Math.abs(xr);
        absx = s * Math.sqrt((xr / s) * (xr / s));
        s = Math.abs(s) + Math.abs(yr);
        absxy = s * Math.sqrt((absx / s) * (absx / s) + (yr / s) * (yr / s));
        P.c = absx / absxy;
        xr = xr / absx;
        P.sr = -(xr * yr) / absxy;
        P.zr = xr * absxy;
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

    public static void genr(double x, double y, PlaneRotationDense P) {

        double s, absx, absxy;

        if (x == 0) {
            P.c = 0.;
            P.sr = 1.;
            P.zr = y;
            return;
        }
        double xr = x;
        double yr = y;

        s = Math.abs(xr);
        absx = s * Math.sqrt((xr / s) * (xr / s));
        s = Math.abs(s) + Math.abs(yr);
        absxy = s * Math.sqrt((absx / s) * (absx / s) + (yr / s) * (yr / s));
        P.c = absx / absxy;
        xr = xr / absx;
        P.sr = -(xr * yr) / absxy;
        P.zr = xr * absxy;
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

    public static PlaneRotationDense genr(MatrixDense A, int ii, int jj1, int jj2) {

        int i = ii;
        int j1 = jj1;
        int j2 = jj2;

        PlaneRotationDense P = PlaneRotationDense.genr(A.unsafeGet(i, j1), A.unsafeGet(i, j2));
        A.unsafeSet(i, j1, P.zr);
        A.unsafeSet(i, j2, 0);
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

    public static void genr(MatrixDense A, int ii, int jj1, int jj2, PlaneRotationDense P) {
        PlaneRotationDense.genr(A.unsafeGet(ii, jj1), A.unsafeGet(ii, jj2), P);
        A.unsafeSet(ii, jj1, P.zr);
        A.unsafeSet(ii, jj2, 0);
        A.unsafeSet(ii, jj2, 0);
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
//    public static PlaneRotationDense genr(double x, double y) {
//
//        PlaneRotationDense P = new PlaneRotationDense();
//
//        double s = Math.abs(x) + Math.abs(y);
//
//        if (s == 0.) {
//            P.c = 1.;
//            P.sr = 0.;
//            P.zr = 0.;
//            return P;
//        }
//
//        P.zr = s * Math.sqrt((x / s) * (x / s) + (y / s) * (y / s));
//        P.c = x / P.zr;
//        P.sr = -y / P.zr;
//        return P;
//    }


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
//    public static void genr(double x, double y, PlaneRotationDense P) {
//
//
//        double s = Math.abs(x) + Math.abs(y);
//
//        if (s == 0.) {
//            P.c = 1.;
//            P.sr = 0.;
//            P.zr = 0.;
//            return;
//        }
//
//        P.zr = s * Math.sqrt((x / s) * (x / s) + (y / s) * (y / s));
//        P.c = x / P.zr;
//        P.sr = -y / P.zr;
//    }


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


    public static void pa(PlaneRotationDense P, MatrixDense A, int ii1, int ii2, int jj1, int jj2) {

        double t1r, t2r;

        for (int j = jj1; j <= jj2; j++) {
            t1r = P.c * A.unsafeGet(ii1, j) + P.sr * A.unsafeGet(ii2, j);
            t2r = P.c * A.unsafeGet(ii2, j) - P.sr * A.unsafeGet(ii1, j);
            A.unsafeSet(ii1, j, t1r);
            A.unsafeSet(ii2, j, t2r);
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
    public static void pha(PlaneRotationDense P, MatrixDense A, int ii1, int ii2, int jj1, int jj2) {

        double t1r, t2r;

        int i1 = ii1;
        int i2 = ii2;
        int j1 = jj1;
        int j2 = jj2;


        for (int j = j1; j <= j2; j++) {
            t1r = P.c * A.unsafeGet(i1, j) - P.sr * A.unsafeGet(i2, j);
            t2r = P.c * A.unsafeGet(i2, j) + P.sr * A.unsafeGet(i1, j);
            A.unsafeSet(i1, j, t1r);
            A.unsafeSet(i2, j, t2r);
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

    public static void ap(MatrixDense A, PlaneRotationDense P, int ii1, int ii2, int jj1, int jj2) {

        double t1r, t1i, t2r, t2i;
        for (int i = ii1; i <= ii2; i++) {
            t1r = P.c * A.unsafeGet(i, jj1) - P.sr * A.unsafeGet(i, jj2);
            t2r = P.c * A.unsafeGet(i, jj2) + P.sr * A.unsafeGet(i, jj1);
            A.unsafeSet(i, jj1, t1r);
            A.unsafeSet(i, jj2, t2r);
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

    public static void aph(MatrixDense A, PlaneRotationDense P, int ii1, int ii2, int jj1, int jj2) {

        double t1r, t1i, t2r, t2i;
        for (int i = ii1; i <= ii2; i++) {
            t1r = P.c * A.unsafeGet(i, jj1) + P.sr * A.unsafeGet(i, jj2);
            t2r = P.c * A.unsafeGet(i, jj2) - P.sr * A.unsafeGet(i, jj1);
            A.unsafeSet(i, jj1, t1r);
            A.unsafeSet(i, jj2, t2r);
        }
    }
}
