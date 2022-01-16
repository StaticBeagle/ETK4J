package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;

import java.util.Arrays;

public class QuadraticSplineZeroCentered extends Spline {

    private static final double P5 = 0.5, P33 = 1.0 / 3.0;

    protected QuadraticSplineZeroCentered(double[] x, double[] y) {
        super(x, 3);
        //TODO check if x is monotonic
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
        if (a < _x[0] || b > _x[_x.length - 1]) {
            throw new IllegalArgumentException(
                    String.format("The spline is not defined outside of [%.4g, %.4g]", _x[0], _x[_x.length - 1]));
        }
        // TODO check that b > a; for all cases look at spline
        int ia = this.findSegmentIndex(a);
        int ib = this.findSegmentIndex(b);

        double integral = evaluateAntiDerivativeAt(ia, _x[ia + 1]) - evaluateAntiDerivativeAt(ia, a);
        for(int i = ia + 1; i < ib; ++i) {
            integral += evaluateAntiDerivativeAt(i, _x[i + 1]) - evaluateAntiDerivativeAt(i, _x[i]);
        }
        return evaluateAntiDerivativeAt(ib, b) - evaluateAntiDerivativeAt(ib, _x[ib]) + integral;
    }

    @Override
    public double integrate(double x) {
        return integrate(_x[0], x);
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
        for (int i = 0, j = 0; i < _x.length - 1; ++i, ++j) {
            double a = coefficients[j];
            double b = coefficients[++j];
            double c = coefficients[++j];

            sb.append("S").append(i + 1).append("(x) = ")
                    .append(a != 0d ? String.format("%.4g * x^2", a, _x[i]) : "")
                    .append(b != 0d ? String.format(" + %.4g * x", b, _x[i]) : "")
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

//        double a = 0;
//        double b = 4;
//        double xi = 3.5;
//
//        LinearSpline ls = LinearSpline.newLinearSpline(x, y);
//        System.out.println(ls.evaluateAt(1.5));
//        System.out.println(ls.differentiate(2.5));
//        System.out.println(ls.integrate(a, b));
//        System.out.println(ls.integrate(xi));
//        System.out.println(ls);
//
//        QuadraticSplineZeroCentered qs = newQuadraticSpline(x, y);
//        System.out.println(qs.evaluateAt(1.5));
//        System.out.println(qs.differentiate(2.5));
//        System.out.println(qs.integrate(a, b));
//        System.out.println(qs.integrate(xi));
//        System.out.println(qs);
//
//        CubicSpline cs = CubicSpline.newNotAKnotSpline(x, y);
//        System.out.println(cs.evaluateAt(1.5));
//        System.out.println(cs.differentiate(2.5));
//        System.out.println(cs.integrate(a, b));
//        System.out.println(cs.integrate(xi));
//        System.out.println(cs);

        x = new double[]{0.0, 10.0, 15.0, 20.0, 22.5, 30.0};
        y = new double[]{0.0, 227.04, 362.78, 517.35, 602.97, 901.67};

//        LinearSpline ls = LinearSpline.newLinearSpline(x, y);
//        System.out.println(ls.evaluateAt(16));
//        System.out.println(ls.differentiate(16));
//        System.out.println(ls.integrate(11, 16));
//        System.out.println(ls);

        QuadraticSplineZeroCentered qs = newQuadraticSpline(x, y);
        System.out.println(qs.evaluateAt(16));
        System.out.println(qs.differentiate(16));
        System.out.println(qs.integrate(11, 16));
        System.out.println(qs.integrate(16));
        System.out.println(qs);

        CubicSpline cs = CubicSpline.newNotAKnotSpline(x, y);
        System.out.println(cs.evaluateAt(16));
        System.out.println(cs.differentiate(16));
        System.out.println(cs.integrate(11, 16));
        System.out.println(cs.integrate(16));
        System.out.println(cs);
    }
}
