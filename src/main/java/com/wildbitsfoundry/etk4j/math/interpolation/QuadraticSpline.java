package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.util.DoubleArrays;

import java.util.Arrays;

public class QuadraticSpline extends Spline {

    private static final double P5 = 0.5, P33 = 1.0 / 3.0;

    protected QuadraticSpline(double[] x, double[] y, double[] dydx) {
        super(x, 3);

        final int n = this.x.length;
        // compute coefficients
        coefficients = new double[(n - 1) * 3]; // 3 coefficients and n - 1 segments
        for (int i = 0, j = 0; i < n - 1; ++i, ++j) {
            double hx = this.x[i + 1] - this.x[i];
            double a = 0.5 * (dydx[i + 1] - dydx[i]) / hx;
            double b = dydx[i];
            double c = y[i];
            coefficients[j] = a;
            coefficients[++j] = b;
            coefficients[++j] = c;
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
    public double evaluateAt(int i, double x) {
        double t = x - this.x[i];
        i *= 3;
        return coefficients[i + 2] + t * (coefficients[i + 1] + t * coefficients[i]);
    }

    @Override
    public double evaluateDerivativeAt(int i, double x) {
        double t = x - this.x[i];
        i *= 3;
        return coefficients[i + 1] + t * 2.0 * coefficients[i];
    }

    @Override
    public double evaluateAntiDerivativeAt(int i, double x) {
        double t = x - this.x[i];
        i *= 3;
        return t * (coefficients[i + 2] + t * (coefficients[i + 1] * P5 + t * coefficients[i] * P33));
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0, j = 0; i < x.length - 1; ++i, ++j) {
            double a = coefficients[j];
            double b = coefficients[++j];
            double c = coefficients[++j];

            sb.append("S").append(i + 1).append("(x) = ")
                    .append(a != 0d ? String.format("%.4g * (x - %.4f)^2", a, x[i]) : "")
                    .append(b != 0d ? String.format(" + %.4g * (x - %.4f)", b, x[i]) : "")
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
//
//        double[] yi = new double[31];
//        for(int i = 0; i <= 30; ++i) {
//            yi[i] = qs2.evaluateAt(i);
//        }

        CubicSpline cs = CubicSpline.newNotAKnotSpline(x, y);

        System.out.println(Arrays.toString(cs.evaluateAt(DoubleArrays.linSteps(0, 30))));
    }
}
