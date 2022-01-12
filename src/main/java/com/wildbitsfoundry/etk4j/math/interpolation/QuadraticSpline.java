package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;
import com.wildbitsfoundry.etk4j.util.NumArrays;

import java.lang.reflect.Array;
import java.util.Arrays;

public class QuadraticSpline extends Spline {

    private static final double P5 = 0.5, P33 = 1.0 / 3.0;

    protected QuadraticSpline(double[] x, double[] y, double[] dydx) {
        super(x, 3);

        final int n = _x.length;
        // compute coefficients
        _coefs = new double[(n - 1) * 3]; // 3 coefficients and n - 1 segments
        for (int i = 0, j = 0; i < n - 1; ++i, ++j) {
            double hx = _x[i + 1] - _x[i];
            if (hx <= 0.0) {
                throw new IllegalArgumentException("x must be monotonically increasing");
            }
            double a = 0.5 * (dydx[i + 1] - dydx[i]) / hx;
            double b = dydx[i];
            double c = y[i];
            _coefs[j] = a;
            _coefs[++j] = b;
            _coefs[++j] = c;
        }
    }

    public static QuadraticSpline newNaturalSpline(double[] x, double[] y) {
        return newNaturalSplineInPlace(Arrays.copyOf(x, x.length), y);
    }

    public static QuadraticSpline newNaturalSplineInPlace(double[] x, double[] y) {
        final int n = x.length;
        double[] b = new double[n];
        // Natural conditions
        b[0] = 0.0;
        for (int i = 1; i < n; ++i) {
            double hx = x[i] - x[i - 1];
            b[i] = 2 * (y[i] - y[i - 1]) / hx - b[i - 1];
        }
        return new QuadraticSpline(x, y, b);
    }

    public static QuadraticSpline newClampedSpline(double[] x, double[] y, double d0) {
        return newClampedSplineInPlace(Arrays.copyOf(x, x.length), y, d0);
    }

    public static QuadraticSpline newClampedSplineInPlace(double[] x, double[] y, double d0) {
        final int n = x.length;
        double[] b = new double[n];
        // Clamped conditions
        b[0] = (y[1] - y[0]) / (x[1] - x[0]) - (x[1] - x[0]) * d0;
        for (int i = 1; i < n; ++i) {
            double hx = x[i] - x[i - 1];
            b[i] = 2 * (y[i] - y[i - 1]) / hx - b[i - 1];
        }
        return new QuadraticSpline(x, y, b);
    }

    @Override
    public double evaluateSegmentAt(int i, double x) {
        double t = x - _x[i];
        i *= 3;
        return _coefs[i + 2] + t * (_coefs[i + 1] + t * _coefs[i]);
    }

    // TODO change t to x
    @Override
    public double evaluateDerivativeAt(int i, double t) {
        i *= 3;
        return _coefs[i + 1] + t * 2.0 * _coefs[i];
    }

    // TODO change t to x
    @Override
    public double evaluateAntiDerivativeAt(int i, double t) {
        i *= 3;
        return t * (_coefs[i + 2] + t * (_coefs[i + 1] * P5 + t * _coefs[i] * P33));
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0, j = 0; i < _x.length - 1; ++i, ++j) {
            double a = _coefs[j];
            double b = _coefs[++j];
            double c = _coefs[++j];

            sb.append("S").append(i + 1).append("(x) = ")
                    .append(a != 0d ? String.format("%.4g * (x - %.4f)^2", a, _x[i]) : "")
                    .append(b != 0d ? String.format(" + %.4g * (x - %.4f)", b, _x[i]) : "")
                    .append(c != 0d ? String.format(" + %.4g", c) : "")
                    .append(System.lineSeparator());
        }
        sb.setLength(Math.max(sb.length() - System.lineSeparator().length(), 0));
        return sb.toString().replace("+ -", "- ").replace("- -", "+ ")
                .replace("=  + ", "= ").replace("=  - ", "= -");
    }

    public static void main(String[] args) {
        double[] x = {0, 1, 2, 3, 4, 5};
        double[] y = {0, 1, 4, 9, 16, 25};


        x = new double[]{0.0, 10.0, 15.0, 20.0, 22.5, 30.0};
        y = new double[]{0.0, 227.04, 362.78, 517.35, 602.97, 901.67};

        QuadraticSpline qs2 = newNaturalSpline(x, y);
        System.out.println(qs2.evaluateAt(16));
        System.out.println(qs2.differentiate(16));
        System.out.println(qs2.integrate(11, 16));
        System.out.println(qs2);

        QuadraticSpline qs3 = newClampedSpline(x, y, 1.0);
        System.out.println(qs3.evaluateAt(16));
        System.out.println(qs3.differentiate(16));
        System.out.println(qs3.integrate(11, 16));
        System.out.println(qs3);

        double[] yi = new double[31];
        for(int i = 0; i <= 30; ++i) {
            yi[i] = qs2.evaluateAt(i);
        }

        System.out.println(Arrays.toString(yi));
    }
}
