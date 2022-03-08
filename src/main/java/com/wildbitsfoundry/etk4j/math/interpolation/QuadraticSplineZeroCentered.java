package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;

import java.util.Arrays;

public class QuadraticSplineZeroCentered extends Spline {

    private static final double P5 = 0.5, P33 = 1.0 / 3.0;

    protected QuadraticSplineZeroCentered(double[] x, double[] y) {
        super(x, 3);
        int dim = x.length - 1;
        Matrix X = new Matrix(3 * dim, 3 * dim);
        for (int i = 0; i < dim; ++i) {
            X.set(2 * (i + 1) - 2, 3 * i + 2, 1.0);
            X.set(2 * (i + 1) - 1, 3 * i + 2, 1.0);
            X.set(2 * (i + 1) - 2, 3 * i, x[i] * x[i]);
            X.set(2 * (i + 1) - 2, 3 * i + 1, x[i]);
            X.set(2 * (i + 1) - 1, 3 * i, x[i + 1] * x[i + 1]);
            X.set(2 * (i + 1) - 1, 3 * i + 1, x[i + 1]);
        }
        for (int i = 0; i < dim - 1; ++i) {
            X.set(2 * dim + i, 3 * i + 1, 1.0);
            X.set(2 * dim + i, 3 * i + 4, -1.0);
            X.set(2 * dim + i, 3 * i, 2 * x[i + 1]);
            X.set(2 * dim + i, 3 * i + 3, -2 * x[i + 1]);
        }
        X.set(3 * dim - 1, 0, 1.0);
        Matrix b = new Matrix(3 * dim, 1);

        b.set(0, 0, y[0]);
        for (int i = 1; i < 2 * dim; ++i) {
            b.set(i, 0, y[(i + 1) / 2]);
        }
        b.set(2 * dim - 1, 0, y[dim]);
        coefficients = X.solve(b).getArray();
    }

    public static QuadraticSplineZeroCentered newQuadraticSpline(double[] x, double[] y) {
        return new QuadraticSplineZeroCentered(Arrays.copyOf(x, x.length), y);
    }

    public static QuadraticSplineZeroCentered newQuadraticSplineInPlace(double[] x, double[] y) {
        return new QuadraticSplineZeroCentered(x, y);
    }

    @Override
    public double integrate(double a, double b) {
        if (a < x[0] || b > x[x.length - 1]) {
            throw new IllegalArgumentException(
                    String.format("The spline is not defined outside of [%.4g, %.4g]", x[0], x[x.length - 1]));
        }
        if(b < a) {
            throw new IllegalArgumentException("The upper integration limit b has to be greater or equal to the lower" +
                    " integration limit a.");
        }
        int ia = this.findSegmentIndex(a);
        int ib = this.findSegmentIndex(b);

        double integral = evaluateAntiDerivativeAt(ia, x[ia + 1]) - evaluateAntiDerivativeAt(ia, a);
        for(int i = ia + 1; i < ib; ++i) {
            integral += evaluateAntiDerivativeAt(i, x[i + 1]) - evaluateAntiDerivativeAt(i, x[i]);
        }
        return evaluateAntiDerivativeAt(ib, b) - evaluateAntiDerivativeAt(ib, x[ib]) + integral;
    }

    @Override
    public double integrate(double x) {
        return integrate(this.x[0], x);
    }

    @Override
    public double evaluateAt(int i, double x) {
        i *= 3;
        return coefficients[i + 2] + x * (coefficients[i + 1] + x * coefficients[i]);
    }

    @Override
    public double evaluateDerivativeAt(int i, double x) {
        i *= 3;
        return coefficients[i + 1] + x * 2.0 * coefficients[i];
    }

    @Override
    public double evaluateAntiDerivativeAt(int i, double x) {
        i *= 3;
        return x * (coefficients[i + 2] + x * (coefficients[i + 1] * P5 + x * coefficients[i] * P33));
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0, j = 0; i < x.length - 1; ++i, ++j) {
            double a = coefficients[j];
            double b = coefficients[++j];
            double c = coefficients[++j];

            sb.append("S").append(i + 1).append("(x) = ")
                    .append(a != 0d ? String.format("%.4g * x^2", a, x[i]) : "")
                    .append(b != 0d ? String.format(" + %.4g * x", b, x[i]) : "")
                    .append(c != 0d ? String.format(" + %.4g", c) : "")
                    .append(System.lineSeparator());
        }
        sb.setLength(Math.max(sb.length() - System.lineSeparator().length(), 0));
        return sb.toString().replace("+ -", "- ").replace("- -", "+ ")
                .replace("=  + ", "= ").replace("=  - ", "= -");
    }
}
