package com.wildbitsfoundry.etk4j.math.curvefitting;

import com.wildbitsfoundry.etk4j.math.polynomial.Polynomial;

import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkXYDimensions;

/**
 * The {@code CurveFitting class} provides various methods for curve fitting. Some of this methods use least squares
 * to approximate the data to be fitted. Most of the fitting results are returned in an array of doubles. This differs
 * from interpolation in the sense that no function is created and/or evaluated but only coefficients are determined.
 */
public final class CurveFitting {
    private CurveFitting() {
    }

    /**
     * Constructs a line <code>f(x) = m * x + b</code> defined by two points (x0, y0) and (x1, y1).
     *
     * @param x0 Starting coordinate for x.
     * @param x1 End coordinate for x.
     * @param y0 Starting coordinate for y.
     * @param y1 End coordinate for y.
     * @return [m, b];
     */
    public static double[] line(double x0, double x1, double y0, double y1) {
        double m = (y1 - y0) / (x1 - x0);
        return new double[]{m, y0 - m * x0};
    }

    /**
     * Linear least squares fit. <br>
     * Fits x and y to a function <code>f(x) = m * x + b</code> in a least square sense.
     *
     * @param x Independent variable. This array must be monotonically increasing.
     * @param y Dependent variable.
     * @return [m, b];
     */
    public static double[] linear(double[] x, double[] y) {
        checkXYDimensions(x, y);

        final int n = x.length;
        double xbar = 0.0;
        double ybar = 0.0;
        double ssxx = 0.0;
        double ssxy = 0.0;
        for (int i = 0; i < n; ++i) {
            xbar += x[i];
            ybar += y[i];
            ssxx += x[i] * x[i];
            ssxy += x[i] * y[i];
        }
        xbar /= n;
        ybar /= n;
        ssxx -= n * xbar * xbar;
        ssxy -= n * xbar * ybar;

        double m = ssxy / ssxx;
        double b = ybar - m * xbar;

        return new double[]{m, b};
    }


    /**
     * Constructs a parabola <code>f(x) = a * x^2 + b * x + c</code> defined by points (x0, y0), (x1, y1), (x2, y2).
     *
     * @param x0 Starting coordinate for x.
     * @param x1 Middle coordinate for x.
     * @param x2 End coordinate for x.
     * @param y0 Starting coordinate for y.
     * @param y1 Middle coordinate for y.
     * @param y2 End coordinate for y.
     * @return [a, b, c];
     */
    public static double[] parabola(double x0, double x1, double x2, double y0, double y1, double y2) {
        double denom = (x0 - x1) * (x0 - x2) * (x1 - x2);
        double a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / denom;
        double b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / denom;
        double c = (x1 * x2 * (x1 - x2) * y0 + x2 * x0 * (x2 - x0) * y1 + x0 * x1 * (x0 - x1) * y2) / denom;

        return new double[]{a, b, c};
    }

    /**
     * Polynomial fit.Fits a polynomial to the input data in exact or least square sense depending on the order and
     * the number of values entered.
     * @param x Array of abscissas. This array must be monotonically increasing.
     * @param y Array of ordinates.
     * @param n Order of the requested polynomial
     * @return The coefficients of the result of the polynomial fit.
     */
    public static double[] polynomial(double[] x, double[] y, int n) {
        return Polynomial.polyFit(x, y, n).getCoefficients();
    }

    /**
     * Exponential least squares fit. <br>
     * Fits x and y to a function a * e^(b * x) in a least square sense.
     *
     * @param x Independent variable. This array must be monotonically increasing.
     * @param y Dependent variable.
     * @return [a, b]
     */
    public static double[] exponential(double[] x, double[] y) {
        checkXYDimensions(x, y);

        final int n = x.length;
        double k0 = 0.0, k1 = 0.0, k2 = 0.0, k3 = 0.0, r0 = 0.0, r1 = 0.0;
        for (int i = 0; i < n; i++) {
            k0 += y[i];
            k3 = x[i] * y[i];
            k1 += k3;
            k2 += k3 * x[i];
            k3 = y[i] * Math.log(y[i]);
            r0 += k3;
            r1 += x[i] * k3;
        }
        double det = 1 / (k0 * k2 - k1 * k1);
        double a = Math.exp((k2 * r0 - k1 * r1) * det);
        double b = (k0 * r1 - k1 * r0) * det;
        return new double[]{a, b};
    }

    /**
     * Logarithmic least squares fit.
     * Fits x and y to a function a + b * ln(x) in a least square sense.
     *
     * @param x Independent variable.
     * @param y Dependent variable.
     * @return [a, b];
     */
    public static double[] logarithmic(double[] x, double[] y) {
        checkXYDimensions(x, y);

        final int n = x.length;
        double k0 = 0.0, k1 = 0.0, k2 = 0.0, k3 = 0.0, k4 = 0.0;
        for (int i = 0; i < n; i++) {
            k0 += y[i];
            k3 = Math.log(x[i]);
            k1 += k3;
            k2 += k3 * k3;
            k4 += y[i] * k3;
        }
        double b = (n * k4 - k0 * k1) / (n * k2 - k1 * k1);
        double a = (k0 - b * k1) / n;
        return new double[]{a, b};
    }

    /**
     * Power Law least squares fit.
     * Fits x and y to a function a * x^(b) in a least square sense.
     *
     * @param x Independent variable
     * @param y Dependent variable
     * @return [a, b];
     */
    public static double[] power(double[] x, double[] y) {
        checkXYDimensions(x, y);

        final int n = x.length;
        double k0 = 0.0, k1 = 0.0, k2 = 0.0, k3 = 0.0, lx, ly;
        for (int i = 0; i < n; i++) {
            lx = Math.log(x[i]);
            ly = Math.log(y[i]);
            k0 += lx;
            k1 += ly;
            k2 += lx * ly;
            k3 += lx * lx;
        }
        double b = (n * k2 - k0 * k1) / (n * k3 - k0 * k0);
        double a = Math.exp((k1 - b * k0) / n);
        return new double[]{a, b};
    }
}
