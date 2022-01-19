package com.wildbitsfoundry.etk4j.math.interpolation;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;

import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkMinXLength;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkXYDimensions;

public class CubicSpline extends Spline {

    private static final double P5 = 0.5, P33 = 1.0 / 3.0, P25 = 0.25;

    // TODO move this to matrix?
    private static void solveLDUTridiagonalSystem(double[] lower, double[] diagonal, double[] upper, double[] b) {
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

    private static double[] solveLDLtTridiagonalSystem(double[] lower, double[] diagonal, double[] b) {
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

    private CubicSpline(double[] x, double[] y, double[] coefficients) {
        super(x, 4);
        this.coefficients = coefficients;
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
        checkXYDimensions(x, y);
        checkMinXLength(x, 3);
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
        solveLDLtTridiagonalSystem(lower, diagonal, r);

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
        int j = coefficients.length - 4;
        coefficients[j] = -r[n - 1] / (6.0 * (x[n + 1] - x[n]));
        coefficients[++j] = r[n - 1] / 2.0;
        coefficients[++j] = (y[n + 1] - y[n]) / (x[n + 1] - x[n]) - (x[n + 1] - x[n]) * (2.0 * r[n - 1]) / 6.0; ;
        coefficients[++j] = y[n];
        return new CubicSpline(x, y, coefficients);
    }

    public static CubicSpline newParabolicallyTerminatedSpline(double[] x, double[] y) {
        return newParabolicallyTerminatedSplineInPlace(Arrays.copyOf(x, x.length), y);
    }

    public static CubicSpline newParabolicallyTerminatedSplineInPlace(double[] x, double[] y) {
        checkXYDimensions(x, y);
        checkMinXLength(x, 3);
        final int n = x.length;
        double[] r = new double[n];

        double[] diagonal = new double[n];
        double[] lower = new double[n - 1];
        double[] upper = new double[n - 1];

        double hn = x[n - 1] - x[n - 2];
        diagonal[0] = 1.0;
        diagonal[n - 1] = -1.0;

        upper[0] = -1.0;
        for (int i = 0; i < n - 2; ++i) {
            double hi = x[i + 1] - x[i];
            double hiPlus1 = x[i + 2] - x[i + 1];
            lower[i] = hi;
            upper[i + 1] = hiPlus1;
            diagonal[i + 1] = 2.0 * (hi + hiPlus1);
            r[i + 1] = 3.0 / hiPlus1 * (y[i + 2] - y[i + 1]) - 3.0 / hi * (y[i + 1] - y[i]);
        }
        lower[lower.length - 1] = 1.0;

        // result is in r
        solveLDUTridiagonalSystem(lower, diagonal, upper, r);

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

    public static CubicSpline newClampedSpline(double[] x, double[] y, double d0, double dn) {
        return newClampedSplineInPlace(Arrays.copyOf(x, x.length), y, d0, dn);
    }

    public static CubicSpline newClampedSplineInPlace(double[] x, double[] y, double d0, double dn) {
        checkXYDimensions(x, y);
        checkMinXLength(x, 2);
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
            r[i + 1] = 3.0 * ((y[i + 2] - y[i + 1]) / hiPlus1  - (y[i + 1] - y[i]) / hi);
        }
        lower[lower.length - 1] = hn;

        // result is in r
        solveLDLtTridiagonalSystem(lower, diagonal, r);

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
        // compute coefficients
        double[] coefficients = new double[(n - 1) * 4]; // 4 coefficients and n - 1 segments
        for (int i = 0, j = 0; i < n - 1; ++i, ++j) {
            double hx = x[i + 1] - x[i];
            double m0 = d[i] * hx;
            double m1 = d[i + 1] * hx;
            coefficients[j] = (2 * y[i] + m0 - 2 * y[i + 1] + m1) / (hx * hx * hx);;
            coefficients[++j] = (-3 * y[i] - 2 * m0 + 3 * y[i + 1] - m1) / (hx * hx);;
            coefficients[++j] = d[i];
            coefficients[++j] = y[i];;
        }
        return new CubicSpline(x, y, coefficients);
    }

    public static CubicSpline newNotAKnotSplineInPlace(double[] x, double[] y) {
        checkXYDimensions(x, y);
        checkMinXLength(x, 4);
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
        solveLDUTridiagonalSystem(lower, diagonal, upper, r);

        double b0 = (1 + h0 / h1) * r[0] - h0 / h1 * r[1];
        double[] coefficients = new double[(x.length - 1) * 4]; // 4 coefficients and n - 1 segments
        coefficients[0] = (r[0] - b0) / (3.0 * (x[1] - x[0]));
        coefficients[1] = b0;
        coefficients[2] = (y[1] - y[0]) / (x[1] - x[0]) - (x[1] - x[0]) * (r[0] + 2 * b0) / 3.0;
        coefficients[3] = y[0];

        final int noPoints = y.length;;
        for (int i = 0, j = 4; i < noPoints - 3; ++i, ++j) {
            coefficients[j] = (r[i + 1] - r[i]) / (3.0 * (x[i + 2] - x[i + 1]));
            coefficients[++j] = r[i];
            coefficients[++j] = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]) - (x[i + 2] - x[i + 1])
                    * (r[i + 1] + 2 * r[i]) / 3.0;
            coefficients[++j] = y[i + 1];
        }

        double bn = (1 + hn / hnMinus1) * r[r.length - 1] - hn / hnMinus1 * r[r.length - 2];
        int j = coefficients.length - 4;
        int m = r.length;
        coefficients[j] = (bn - r[m - 1]) / (3.0 * (x[noPoints - 1] - x[noPoints - 2]));
        coefficients[++j] = r[m - 1];
        coefficients[++j] = (y[noPoints - 1] - y[noPoints - 2]) / (x[noPoints - 1] - x[noPoints - 2]) -
                (x[noPoints - 1] - x[noPoints - 2]) * (bn + 2 * r[m - 1]) / 3.0;
        coefficients[++j] = y[noPoints - 2];
        return new CubicSpline(x, y, coefficients);
    }

    public static CubicSpline newNotAKnotSpline(double[] x, double[] y) {
        return newNotAKnotSplineInPlace(Arrays.copyOf(x, x.length), y);
    }

    @Override
    public double evaluateAt(int i, double x) {
        double t = x - this.x[i];
        i <<= 2;
        return coefficients[i + 3] + t * (coefficients[i + 2] + t * (coefficients[i + 1] + t * coefficients[i]));
    }

    @Override
    public double evaluateDerivativeAt(int i, double x) {
        double t = x - this.x[i];
        i <<= 2;
        return coefficients[i + 2] + t * (2 * coefficients[i + 1] + t * 3 * coefficients[i]);
    }

    @Override
    public double evaluateAntiDerivativeAt(int i, double x) {
        double t = x - this.x[i];
        i <<= 2;
        return t * (coefficients[i + 3] + t * (coefficients[i + 2] * P5 + t * (coefficients[i + 1] * P33 + t * coefficients[i] * P25)));
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0, j = 0; i < x.length - 1; ++i, ++j) {
            double a = coefficients[j];
            double b = coefficients[++j];
            double c = coefficients[++j];
            double d = coefficients[++j];

            sb.append("S").append(i + 1).append("(x) = ")
                    .append(a != 0d ? String.format("%.4g * (x - %.4f)^3", a, x[i]) : "")
                    .append(b != 0d ? String.format(" + %.4g * (x - %.4f)^2", b, x[i]) : "")
                    .append(c != 0d ? String.format(" + %.4g * (x - %.4f)", c, x[i]) : "")
                    .append(d != 0d ? String.format(" + %.4g", d, x[i]) : "").append(System.lineSeparator());
        }
        sb.setLength(Math.max(sb.length() - System.lineSeparator().length(), 0));
        return sb.toString().replace("+ -", "- ").replace("- -", "+ ")
                .replace("=  + ", "= ").replace("=  - ", "= -");
    }
}