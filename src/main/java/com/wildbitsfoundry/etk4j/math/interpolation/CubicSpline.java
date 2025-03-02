package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.constant.ConstantsETK;

import java.util.Arrays;

import static com.wildbitsfoundry.etk4j.math.linearalgebra.TridiagonalSolver.solveLDLtTridiagonalSystem;
import static com.wildbitsfoundry.etk4j.math.linearalgebra.TridiagonalSolver.solveLDUTridiagonalSystem;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkMinXLength;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkXYDimensions;

public class CubicSpline extends Spline {

    private static final double P5 = 0.5, P33 = 1.0 / 3.0, P25 = 0.25;

    private CubicSpline(double[] x, double[] coefficients) {
        super(x, 4);
        this.coefficients = coefficients;
    }

    /**
     * Creates a new {@code CubicSpline} with not-a-knot conditions. This method is an alias for
     * {@link #newNotAKnotSpline(double[], double[])}.
     *
     * @param x The array of x coordinates. The values in this array must be unique and strictly increasing.
     *          A copy of this array is made internally. This array must contain 4 values or more.
     * @param y The array of y coordinates.
     * @return A new instance of a cubic spline with not-a-knot conditions.
     */
    public static CubicSpline newCubicSpline(double[] x, double[] y) {
        return newNotAKnotSpline(x, y);
    }

    /**
     * Creates a new {@code CubicSpline} with not-a-knot conditions. This method is an alias for
     * {@link #newNotAKnotSplineInPlace(double[], double[])}.
     *
     * @param x The array of x coordinates. The values in this array must be unique and strictly increasing.This array
     *          is not copied so any changes to this array will be reflected in the spline. This array must contain 4
     *          values or more.
     * @param y The array of y coordinates.
     * @return A new instance of a cubic spline with not-a-knot conditions.
     */
    public static CubicSpline newCubicSplineInPlace(double[] x, double[] y) {
        return newNotAKnotSplineInPlace(x, y);
    }

    /**
     * Creates a new {@code CubicSpline} with natural conditions.
     *
     * @param x The array of x coordinates. The values in this array must be unique and strictly increasing.
     *          A copy of this array is made internally. This array must contain 3 values or more.
     * @param y The array of y coordinates.
     * @return A new instance of a cubic spline with natural conditions.
     * @see <a href="https://www.bragitoff.com/2018/02/cubic-spline-piecewise-interpolation-c-program/">Natural Spline</a>
     */
    public static CubicSpline newNaturalSpline(double[] x, double[] y) {
        return newNaturalSplineInPlace(Arrays.copyOf(x, x.length), y);
    }

    /**
     * Creates a new {@code CubicSpline} with natural conditions.
     *
     * @param x The array of x coordinates. The values in this array must be unique and strictly increasing.This array
     *          is not copied so any changes to this array will be reflected in the spline. This array must contain 3
     *          values or more.
     * @param y The array of y coordinates.
     * @return A new instance of a cubic spline with natural conditions.
     * @see <a href="https://www.bragitoff.com/2018/02/cubic-spline-piecewise-interpolation-c-program/">Natural Spline</a>
     */
    public static CubicSpline newNaturalSplineInPlace(double[] x, double[] y) {
        checkXYDimensions(x, y);
        checkMinXLength(x, 3);
        final int n = x.length - 2;
        double[] diagonal = new double[n];
        double[] lower = new double[n - 1];
        double[] r = new double[n];
        for (int i = 0; i < n - 1; ++i) {
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
        for (int i = 1, j = 4; i < x.length - 2; ++i, ++j) {
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
        coefficients[++j] = (y[n + 1] - y[n]) / (x[n + 1] - x[n]) - (x[n + 1] - x[n]) * (2.0 * r[n - 1]) / 6.0;
        coefficients[++j] = y[n];
        return new CubicSpline(x, coefficients);
    }

    /**
     * Creates a new {@code CubicSpline} with parabolically terminated conditions.
     *
     * @param x The array of x coordinates. The values in this array must be unique and strictly increasing.
     *          A copy of this array is made internally. This array must contain 3 values or more.
     * @param y The array of y coordinates.
     * @return A new instance of a cubic spline with parabolically terminated conditions.
     */
    public static CubicSpline newParabolicallyTerminatedSpline(double[] x, double[] y) {
        return newParabolicallyTerminatedSplineInPlace(Arrays.copyOf(x, x.length), y);
    }

    /**
     * Creates a new {@code CubicSpline} with parabolically terminated conditions.
     *
     * @param x The array of x coordinates. The values in this array must be unique and strictly increasing.This array
     *          is not copied so any changes to this array will be reflected in the spline. This array must contain 3
     *          values or more.
     * @param y The array of y coordinates.
     * @return A new instance of a cubic spline with parabolically terminated conditions.
     */
    public static CubicSpline newParabolicallyTerminatedSplineInPlace(double[] x, double[] y) {
        checkXYDimensions(x, y);
        checkMinXLength(x, 3);
        final int n = x.length;
        double[] r = new double[n];

        double[] diagonal = new double[n];
        double[] lower = new double[n - 1];
        double[] upper = new double[n - 1];

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
        for (int i = 0, j = 0; i < n - 1; ++i, ++j) {
            double hi = x[i + 1] - x[i];
            double dy = y[i + 1] - y[i];
            coefficients[j] = (r[i + 1] - r[i]) / (3.0 * hi);
            coefficients[++j] = r[i];
            coefficients[++j] = dy / hi - hi * (r[i + 1] + 2.0 * r[i]) / 3.0;
            coefficients[++j] = y[i];
        }
        return new CubicSpline(x, coefficients);
    }

    /**
     * Creates a new {@code CubicSpline} with clamped conditions.
     *
     * @param x  The array of x coordinates. The values in this array must be unique and strictly increasing.
     *           A copy of this array is made internally. This array must contain 2 values or more.
     * @param y  The array of y coordinates.
     * @param d0 The value of the derivative at the first endpoint.
     * @param dn The value of the derivative at the last endpoint.
     * @return A new instance of a cubic spline with clamped conditions.
     */
    public static CubicSpline newClampedSpline(double[] x, double[] y, double d0, double dn) {
        return newClampedSplineInPlace(Arrays.copyOf(x, x.length), y, d0, dn);
    }

    /**
     * Creates a new {@code CubicSpline} with clamped conditions.
     *
     * @param x  The array of x coordinates. The values in this array must be unique and strictly increasing.This array
     *           is not copied so any changes to this array will be reflected in the spline. This array must contain 2
     *           values or more.
     * @param y  The array of y coordinates.
     * @param d0 The value of the derivative at the first endpoint.
     * @param dn The value of the derivative at the last endpoint.
     * @return A new instance of a cubic spline with clamped conditions.
     */
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
            r[i + 1] = 3.0 * ((y[i + 2] - y[i + 1]) / hiPlus1 - (y[i + 1] - y[i]) / hi);
        }
        lower[lower.length - 1] = hn;

        // result is in r
        solveLDLtTridiagonalSystem(lower, diagonal, r);

        double[] coefficients = new double[(x.length - 1) * 4]; // 4 coefficients and n - 1 segments
        for (int i = 0, j = 0; i < n - 1; ++i, ++j) {
            double hi = x[i + 1] - x[i];
            double dy = y[i + 1] - y[i];
            coefficients[j] = (r[i + 1] - r[i]) / (3.0 * hi);
            coefficients[++j] = r[i];
            coefficients[++j] = dy / hi - hi * (r[i + 1] + 2.0 * r[i]) / 3.0;
            coefficients[++j] = y[i];
        }
        return new CubicSpline(x, coefficients);
    }

    // <copyright file="CubicSpline.cs" company="Math.NET">
    // Math.NET Numerics, part of the Math.NET Project
    // http://numerics.mathdotnet.com
    // http://github.com/mathnet/mathnet-numerics
    // See https://github.com/StaticBeagle/ETK4J/blob/master/Math.NET.txt.

    /**
     * Creates a new Akima {@code CubicSpline}.
     *
     * @param x The array of x coordinates. The values in this array must be unique and strictly increasing.
     *          A copy of this array is made internally. This array must contain 5 values or more.
     * @param y The array of y coordinates.
     * @return A new instance of an Akima cubic spline.
     */
    public static CubicSpline newAkimaSpline(double[] x, double[] y) {
        return newAkimaSplineInPlace(x, y);
    }

    /**
     * Creates a new Akima {@code CubicSpline}.
     *
     * @param x The array of x coordinates. The values in this array must be unique and strictly increasing.This array
     *          is not copied so any changes to this array will be reflected in the spline. This array must contain 5
     *          values or more.
     * @param y The array of y coordinates.
     * @return A new instance of an Akima cubic spline.
     */
    public static CubicSpline newAkimaSplineInPlace(double[] x, double[] y) {
        checkXYDimensions(x, y);
        checkMinXLength(x, 5);
        /* Prepare divided differences (diff) and weights (w) */

        double[] diff = new double[x.length - 1];
        double[] weights = new double[x.length - 1];

        for (int i = 0; i < diff.length; i++) {
            diff[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
        }

        for (int i = 1; i < weights.length; i++) {
            weights[i] = Math.abs(diff[i] - diff[i - 1]);
        }

        /* Prepare Hermite interpolation scheme */

        double[] dd = new double[x.length];

        for (int i = 2; i < dd.length - 2; i++) {
            dd[i] = weights[i - 1] < ConstantsETK.DOUBLE_EPS && weights[i + 1] < ConstantsETK.DOUBLE_EPS
                    ? (((x[i + 1] - x[i]) * diff[i - 1]) + ((x[i] - x[i - 1]) * diff[i])) / (x[i + 1] - x[i - 1])
                    : ((weights[i + 1] * diff[i - 1]) + (weights[i - 1] * diff[i])) / (weights[i + 1] + weights[i - 1]);
        }

        dd[0] = differentiateThreePoint(x, y, 0, 0, 1, 2);
        dd[1] = differentiateThreePoint(x, y, 1, 0, 1, 2);
        dd[x.length - 2] = differentiateThreePoint(x, y, x.length - 2, x.length - 3, x.length - 2, x.length - 1);
        dd[x.length - 1] = differentiateThreePoint(x, y, x.length - 1, x.length - 3, x.length - 2, x.length - 1);

        // compute coefficients
        double[] coefficients = new double[(x.length - 1) * 4]; // 4 coefficients and n - 1 segments
        for (int i = 0, j = 0; i < x.length - 1; ++i, ++j) {
            double w = x[i + 1] - x[i];
            double w2 = w * w;
            coefficients[j] = (2 * (y[i] - y[i + 1]) / w + dd[i] + dd[i + 1]) / w2;
            coefficients[++j] = (3 * (y[i + 1] - y[i]) / w - 2 * dd[i] - dd[i + 1]) / w;
            coefficients[++j] = dd[i];
            coefficients[++j] = y[i];
        }
        return new CubicSpline(x, coefficients);
    }

    private static double differentiateThreePoint(double[] xx, double[] yy, int indexT, int index0, int index1, int index2) {
        double x0 = yy[index0];
        double x1 = yy[index1];
        double x2 = yy[index2];

        double t = xx[indexT] - xx[index0];
        double t1 = xx[index1] - xx[index0];
        double t2 = xx[index2] - xx[index0];

        double a = (x2 - x0 - (t2 / t1 * (x1 - x0))) / (t2 * (t2 - t1));
        double b = (x1 - x0 - a * t1 * t1) / t1;
        return (2 * a * t) + b;
    }

    /**
     * Creates a new {@code CubicSpline} with not-a-knot conditions
     *
     * @param x The array of x coordinates. The values in this array must be unique and strictly increasing.
     *          A copy of this array is made internally. This array must contain 4 values or more.
     * @param y The array of y coordinates.
     * @return A new instance of a cubic spline with not-a-knot conditions.
     */
    public static CubicSpline newNotAKnotSpline(double[] x, double[] y) {
        return newNotAKnotSplineInPlace(Arrays.copyOf(x, x.length), y);
    }

    /**
     * Creates a new {@code CubicSpline} with not-a-knot conditions
     *
     * @param x The array of x coordinates. The values in this array must be unique and strictly increasing.This array
     *          is not copied so any changes to this array will be reflected in the spline. This array must contain 4
     *          values or more.
     * @param y The array of y coordinates.
     * @return A new instance of a cubic spline with not-a-knot conditions.
     */
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

        final int noPoints = y.length;
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
        return new CubicSpline(x, coefficients);
    }

    /**
     * Evaluate the spline at a given index.
     *
     * @param index The index of the segment to be evaluated.
     * @param x     The value at which to evaluate the spline.
     * @return The result of evaluating the spline at {@code x}.
     */
    @Override
    public double evaluateAt(int index, double x) {
        double t = x - this.x[index];
        index <<= 2;
        return coefficients[index + 3] + t * (coefficients[index + 2] + t * (coefficients[index + 1] + t * coefficients[index]));
    }

    /**
     * Evaluate the derivative of the spline at a given index.
     *
     * @param index The index of the segment to be evaluated.
     * @param x     The value at which to evaluate the derivative of the spline.
     * @return The result of evaluating the derivative of the spline at {@code x}.
     */
    @Override
    public double evaluateDerivativeAt(int index, double x) {
        double t = x - this.x[index];
        index <<= 2;
        return coefficients[index + 2] + t * (2 * coefficients[index + 1] + t * 3 * coefficients[index]);
    }

    /**
     * Evaluate the antiderivative of the spline at a given index.
     *
     * @param index The index of the segment to be evaluated.
     * @param x     The value at which to evaluate the antiderivative of the spline.
     * @return The result of evaluating the antiderivative of the spline at {@code x}.
     */
    @Override
    public double evaluateAntiDerivativeAt(int index, double x) {
        double t = x - this.x[index];
        index <<= 2;
        return t * (coefficients[index + 3] + t * (coefficients[index + 2] * P5 + t * (coefficients[index + 1] * P33 + t * coefficients[index] * P25)));
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
                    .append(d != 0d ? String.format(" + %.4g", d) : "").append(System.lineSeparator());
        }
        sb.setLength(Math.max(sb.length() - System.lineSeparator().length(), 0));
        return sb.toString().replace("+ -", "- ").replace("- -", "+ ")
                .replace("=  + ", "= ").replace("=  - ", "= -");
    }
}