package com.wildbitsfoundry.etk4j.math.interpolation;

import java.lang.reflect.Array;
import java.util.Arrays;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;
import com.wildbitsfoundry.etk4j.util.NumArrays;

import javax.naming.PartialResultException;

import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkMinXLength;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkXYDimensions;

public class CubicSpline extends Spline {

    private static final double P5 = 0.5, P33 = 1.0 / 3.0, P25 = 0.25;

    // TODO move this to matrix?
    private static void solveTridiagonalSystem(double[] lower, double[] diagonal, double[] upper, double[] b) {
        final int m = b.length;
        b[0] = b[0] / diagonal[0];
        // Forward Substitution
        for (int i = 0; i < m - 1; ++i) {
            upper[i] = upper[i] / diagonal[i];
            diagonal[i + 1] = diagonal[i + 1] - lower[i] * upper[i];
            b[i + 1] = (b[i + 1] - lower[i] * b[i]) / diagonal[i + 1];
        }
        // Backwards Substitution
        for (int i = m - 2; i >= 0; --i) {
            b[i] = b[i] - upper[i] * b[i + 1];
        }
    }

    private static double[] solveLDLTSystem(double[] lower, double[] diagonal, double[] b) {
        final int length = diagonal.length;
        for (int i = 0; i < length - 1; ++i) {
            double ui = lower[i];
            lower[i] /= diagonal[i];
            diagonal[i + 1] -= ui * lower[i];
            b[i + 1] -= lower[i] * b[i];
        }

        b[length - 1] /= diagonal[length - 1];
        for (int i = length - 2; i >= 0; --i) {
            b[i] = b[i] / diagonal[i] - lower[i] * b[i + 1];
        }
        return b;
    }

    // TODO check monotonicity
    private CubicSpline(double[] x, double[] y, double[] coefs) {
        super(x, 4);
        _coefs = coefs;
    }

    public static CubicSpline newCubicSpline(double[] x, double[] y) {
        return newNotAKnotSpline(x, y);
    }

    public static CubicSpline newCubicSplineInPlace(double[] x, double[] y) {
        return newNotAKnotSplineInPlace(x, y);
    }

    public static CubicSpline newNaturalSpline(double[] x, double[] y) {
        return newNaturalSplineInPlace(Arrays.copyOf(x, x.length), y);
    }

    public static CubicSpline newNaturalSplineInPlace(double[] x, double[] y) {
        final int n = x.length - 2;
        double[] diagonal = new double[n];
        double[] lower = new double[n - 1];
        double[] r = new double[n];
        for(int i = 0; i < n - 1; ++i) {
            double hi = x[i + 1] - x[i];
            double hiPlus1 = x[i + 2] - x[i + 1];
            double si = (y[i + 1] - y[i]) / hi;
            double siPlus1 = (y[i + 2] - y[i + 1]) / hiPlus1;

            diagonal[i] = 2.0 * (hi + hiPlus1);
            lower[i] = hiPlus1;
            r[i] = 6.0 * (siPlus1 - si);
        }
        double hn = x[n] - x[n - 1];
        double hiPlus1 = x[n + 1] - x[n];
        double sn = (y[n] - y[n - 1]) / hn;
        double snPlus1 = (y[n + 1] - y[n]) / hiPlus1;

        diagonal[n - 1] = 2.0 * (hn + hiPlus1);
        r[n - 1] = 6.0 * (snPlus1 - sn);
        solveLDLTSystem(lower, diagonal, r);

        double[] coefficients = new double[(x.length - 1) * 4];
        coefficients[0] = r[0] / (6.0 * (x[1] - x[0]));
        coefficients[1] = 0;
        coefficients[2] = (y[1] - y[0]) / (x[1] - x[0]) - (x[1] - x[0]) * r[0] / 6.0;
        coefficients[3] = y[0];
        for(int i = 1, j = 4; i < x.length - 2; ++i, ++j) {
            double hi = x[i + 1] - x[i];
            double dy = y[i + 1] - y[i];
            coefficients[j] = (r[i] - r[i - 1]) / (6.0 * hi);
            coefficients[++j] = r[i - 1] / 2.0;
            coefficients[++j] = dy / hi - hi * (2.0 * r[i - 1] + r[i]) / 6.0;
            coefficients[++j] = y[i];
        }
        int m = coefficients.length - 4;
        coefficients[m] = -r[n - 1] / (6.0 * (x[n + 1] - x[n]));
        coefficients[++m] = r[n - 1] / 2.0;
        coefficients[++m] = (y[n + 1] - y[n]) / (x[n + 1] - x[n]) - (x[n + 1] - x[n]) * (2.0 * r[n - 1]) / 6.0; ;
        coefficients[++m] = y[n];
        return new CubicSpline(x, y, coefficients);
    }

    //TODO
//    public static CubicSpline newParabolicallyTerminatedSpline(double[] x, double[] y) {
//        return newParabolicallyTerminatedSplineInPlace(Arrays.copyOf(x, x.length), y);
//    }
//
//    public static double[] computeParabolicallyTerminatedSplineDerivatives(double[] x, double[] y) {
//        checkMinXLength(x, 2);
//        TridiagonalLDLTSystem T = setupLDLT(x, y);
//        final int n = x.length;
//        T.D[0] = T.L[0];
//        double r0 = (y[1] - y[0]) * T.L[0] * T.L[0];
//        T.b[0] = 2.0 * r0;
//        T.D[n - 1] = T.L[n - 2];
//        r0 = (y[n - 1] - y[n - 2]) * T.L[n - 2] * T.L[n - 2];
//        T.b[n - 1] = 2.0 * r0;
//        return T.solve();
//    }

//    public static CubicSpline newParabolicallyTerminatedSplineInPlace(double[] x, double[] y) {
//        double[] dydx = computeParabolicallyTerminatedSplineDerivatives(x, y);
//        return new CubicSpline(x, y, dydx);
//    }

    // TODO
//    public static double[] computeClampedSplineDerivatives(double[] x, double[] y, double d0, double dn) {
//        checkMinXLength(x, 2);
//        TridiagonalLDLTSystem T = setupLDLT(x, y);
//        final int n = x.length;
//        T.b[0] = d0;
//        T.b[T.b.length - 1] = dn;
//        T.b[1] = T.b[1] - T.b[0] * T.L[0];
//        T.b[n - 2] = T.b[n - 2] - T.b[n - 1] * T.L[n - 2];
//        return T.solve(T.b.length - 1, 1);
//    }
//
    public static CubicSpline newClampedSpline(double[] x, double[] y, double d0, double dn) {
        return newClampedSplineInPlace(Arrays.copyOf(x, x.length), y, d0, dn);
    }

    public static CubicSpline newClampedSplineInPlace(double[] x, double[] y, double d0, double dn) {
        final int n = x.length;
        double[] r = new double[n];

        double[] diagonal = new double[n];
        double[] lower = new double[n - 1];

        double h0 = x[1] - x[0];
        double hn = x[n - 1] - x[n - 2];
        diagonal[0] = 2.0 * h0;
        diagonal[n - 1] = 2.0 * hn;

        r[0] = 3.0 / h0 * (y[1] - y[0]) - 3.0 * d0;
        r[n - 1] = 3.0 * dn - 3.0 / hn * (y[n - 1] - y[n - 2]);
        for (int i = 0; i < n - 2; ++i) {
            double hi = x[i + 1] - x[i];
            double hiPlus1 = x[i + 2] - x[i + 1];
            lower[i] = hi;
            diagonal[i + 1] = 2.0 * (hi + hiPlus1);
            r[i + 1] = 3.0 / hiPlus1 * (y[i + 2] - y[i + 1]) - 3.0 / hi * (y[i + 1] - y[i]);
        }
        lower[lower.length - 1] = hn;

        // result is in r
        solveLDLTSystem(lower, diagonal, r);

        double[] coefficients = new double[(x.length - 1) * 4]; // 4 coefficients and n - 1 segments
        for(int i = 0, j = 0; i < n - 1; ++i, ++j) {
            double hi = x[i + 1] - x[i];
            double dy = y[i + 1] - y[i];
            coefficients[j] = (r[i + 1] - r[i]) / (3.0 * hi);
            coefficients[++j] = r[i];
            coefficients[++j] = dy / hi - hi * (r[i + 1] + 2.0 * r[i]) / 3.0;
            coefficients[++j] = y[i];
        }
        return new CubicSpline(x, y, coefficients);
    }

    public static CubicSpline newAkimaSpline(double[] x, double[] y) {
        return newAkimaSpline(x, y, ConstantsETK.DOUBLE_EPS);
    }

    public static CubicSpline newAkimaSplineInPlace(double[] x, double[] y) {
        return newAkimaSplineInPlace(x, y, ConstantsETK.DOUBLE_EPS);
    }

    public static CubicSpline newAkimaSpline(double[] x, double[] y, double ep) {
        return newAkimaSplineInPlace(Arrays.copyOf(x, x.length), y, ep);
    }

    public static CubicSpline newAkimaSplineInPlace(double[] x, double[] y, double ep) {
        checkXYDimensions(x, y);
        checkMinXLength(x, 5);
        final int n = x.length;
        double[] d = new double[n];
        double[] t = new double[n + 3];

        // inner knots
        for (int i = 0; i < n - 1; ++i) {
            t[i + 2] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
        }
        // end points
        t[n + 1] = 2.0 * t[n] - t[n - 1];
        t[n + 2] = 2.0 * t[n + 1] - t[n];
        t[1] = 2.0 * t[2] - t[3];
        t[0] = 2.0 * t[1] - t[2];

        for (int i = 0; i < n; ++i) {
            double c0 = Math.abs(t[i + 3] - t[i + 2]);
            double c1 = Math.abs(t[i + 1] - t[i]);
            double c2 = c0 + c1;
            if (c2 > ep) {
                d[i] = (c0 * t[i + 1] + c1 * t[i + 2]) / c2;
            } else {
                d[i] = 0.5 * (t[i + 2] + t[i + 1]);
            }
        }
        return new CubicSpline(x, y, d);
    }


    public static CubicSpline newNotAKnotSplineInPlace(double[] x, double[] y) {
        final int n = x.length - 2;
        double[] r = new double[n];

        double[] upper = new double[n - 1];
        double[] diagonal = new double[n];
        double[] lower = new double[n - 1];

        for (int i = 0; i < n - 1; ++i) {
            double hiMinus1 = x[i + 1] - x[i];
            double hi = x[i + 2] - x[i + 1];
            double deltaYMinus1 = y[i + 1] - y[i];
            double deltaYi = y[i + 2] - y[i + 1];

            upper[i] = hi;
            diagonal[i] = 2 * (hiMinus1 + hi);
            lower[i] = hi;
            r[i] = 3 * deltaYi / hi - 3 * deltaYMinus1 / hiMinus1;
        }
        diagonal[n - 1] = 2 * (x[n + 1] - x[n - 1]);
        double h0 = x[1] - x[0];
        double h1 = x[2] - x[1];
        diagonal[0] += h0 + h0 * h0 / h1;

        upper[0] -= h0 * h0 / h1;

        double hn = x[n + 1] - x[n];
        double hnMinus1 = x[n] - x[n - 1];
        diagonal[n - 1] += hn + hn * hn / hnMinus1;

        lower[n - 2] -= hn * hn / hnMinus1;

        double deltaYn = y[n + 1] - y[n];
        double deltaYnMinus1 = y[n] - y[n - 1];
        r[n - 1] = 3.0 * deltaYn / hn - 3.0 * deltaYnMinus1 / hnMinus1;
        // result is in r
        solveTridiagonalSystem(lower, diagonal, upper, r);


        double[] b = new double[r.length + 2];
        System.arraycopy(r, 0, b, 1, r.length);
        b[0] = (1 + h0 / h1) * b[1] - h0 / h1 * b[2];
        b[b.length - 1] = (1 + hn / hnMinus1) * b[b.length - 2] - hn / hnMinus1 * b[b.length - 3];

        double[] coefficients = new double[(x.length - 1) * 4]; // 4 coefficients and n - 1 segments
        for (int i = 0, j = 0; i < y.length - 1; ++i, ++j) {
            coefficients[j] = (b[i + 1] - b[i]) / (3.0 * (x[i + 1] - x[i]));
            coefficients[++j] = b[i];
            coefficients[++j] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (x[i + 1] - x[i]) * (b[i + 1] + 2 * b[i]) / 3.0;
            coefficients[++j] = y[i];
        }
        return new CubicSpline(x, y, coefficients);
    }

    public static CubicSpline newNotAKnotSpline(double[] x, double[] y) {
        return newNotAKnotSplineInPlace(Arrays.copyOf(x, x.length), y);
    }

    @Override
    public double evaluateSegmentAt(int i, double x) {
        double t = x - _x[i];
        i <<= 2;
        return _coefs[i + 3] + t * (_coefs[i + 2] + t * (_coefs[i + 1] + t * _coefs[i]));
    }

    @Override
    public double evaluateDerivativeAt(int i, double t) {
        i <<= 2;
        return _coefs[i + 2] + t * (2 * _coefs[i + 1] + t * 3 * _coefs[i]);
    }

    @Override
    public double evaluateAntiDerivativeAt(int i, double t) {
        i <<= 2;
        return t * (_coefs[i + 3] + t * (_coefs[i + 2] * P5 + t * (_coefs[i + 1] * P33 + t * _coefs[i] * P25)));
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0, j = 0; i < _x.length - 1; ++i, ++j) {
            double a = _coefs[j];
            double b = _coefs[++j];
            double c = _coefs[++j];
            double d = _coefs[++j];

            sb.append("S").append(i + 1).append("(x) = ")
                    .append(a != 0d ? String.format("%.4g * (x - %.4f)^3", a, _x[i]) : "")
                    .append(b != 0d ? String.format(" + %.4g * (x - %.4f)^2", b, _x[i]) : "")
                    .append(c != 0d ? String.format(" + %.4g * (x - %.4f)", c, _x[i]) : "")
                    .append(d != 0d ? String.format(" + %.4g", d, _x[i]) : "").append(System.lineSeparator());
        }
        sb.setLength(Math.max(sb.length() - System.lineSeparator().length(), 0));
        return sb.toString().replace("+ -", "- ").replace("- -", "+ ")
                .replace("=  + ", "= ").replace("=  - ", "= -");
    }

    public static void main(String[] args) {
        double[] x = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        double[] y = {0, 1, 4, 9, 16, 25, 36, 49, 64};

        x = new double[]{1, 3, 8, 14, 15, 20, 26, 36, 45, 62, 95};
        y = new double[]{10, 12, 15, 24, 14, 10, 10.5, 15, 50, 60, 85};

//        x = new double[]{1, 3, 8, 14, 15};
//        y = new double[]{10, 12, 15, 24, 14};

//        x = new double[]{0.9, 1.3, 1.9, 2.1};
//        y = new double[]{1.3, 1.5, 1.85, 2.1};

//        CubicSpline cs = CubicSpline.newNotAKnotSpline(x, y);
//        System.out.println(cs.differentiate(2.0));
//        System.out.println(cs.evaluateAntiDerivativeAt(0, 3.0));
//        System.out.println(cs.integrate(3.0));
//        System.out.println(cs.integrate(1.0, 3.0));
//        System.out.println(cs.evaluateAt(7.5));
//        System.out.println(cs);

        CubicSpline csc = CubicSpline.newNaturalSplineInPlace(x, y);
        System.out.println(csc.differentiate(2.0));
        System.out.println(csc.evaluateAntiDerivativeAt(0, 3.0));
        System.out.println(csc.integrate(3.0));
        System.out.println(csc.integrate(1.0, 3.0));
        System.out.println(csc.evaluateAt(7.5));
        System.out.println(csc);


    }
}