package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;
import com.wildbitsfoundry.etk4j.util.NumArrays;

import java.lang.reflect.Array;
import java.util.Arrays;

public class QuadraticSpline extends Spline {

    private static final double P5 = 0.5, P33 = 1.0 / 3.0;

    private static class TridiagonalSystem {
        public double[] L; // Sub-diagonal
        public double[] D; // Diagonal
        public double[] U; // Super-diagonal
        public double[] b; // Solution vector

        public TridiagonalSystem(int n) {
            L = new double[n - 1];
            D = new double[n];
            U = new double[n - 1];
            b = new double[n];
        }

        public double[] solve() {
            int n = b.length;
            for (int i = 0; i < n - 1; ++i) {
                L[i] /= D[i];
                D[i + 1] -= U[i] * L[i];
                b[i + 1] -= L[i] * b[i];
            }

            b[n - 1] /= D[n - 1];
            for (int i = n - 2; i >= 0; --i) {
                b[i] = (b[i] - U[i] * b[i + 1]) / D[i];
            }
            return b;
        }
    }

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

        TridiagonalSystem T = setupSpline(x, y);
        // Natural conditions
        T.b[0] = 0.0;
        T.D[0] = 1.0;

        return new QuadraticSpline(x, y, T.solve());
    }

    public static QuadraticSpline newClampedSpline(double[] x, double[] y, double d0) {
        return newClampedSplineInPlace(Arrays.copyOf(x, x.length), y, d0);
    }

    public static QuadraticSpline newClampedSplineInPlace(double[] x, double[] y, double d0) {

        TridiagonalSystem T = setupSpline(x, y);
        // Natural conditions
        T.b[0] = d0;
        T.D[0] = 1.0;

        return new QuadraticSpline(x, y, T.solve());
    }

    private static TridiagonalSystem setupSpline(double[] x, double[] y) {
        final int n = x.length;
        TridiagonalSystem T = new TridiagonalSystem(n);

        // U is always zero
        // L and D are always one
        Arrays.fill(T.L, 1.0);
        Arrays.fill(T.D, 1.0);

        for (int i = 1; i < n; ++i) {
            double hx = x[i] - x[i - 1];
            T.b[i] = 2 * (y[i] - y[i - 1]) / hx;
        }
        return T;
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
//
//        LinearSpline ls = LinearSpline.newLinearSpline(x, y);
//        System.out.println(ls.evaluateAt(2.5));
//        System.out.println(ls.differentiate(2.5));
//        System.out.println(ls.integrate(3));
//        System.out.println(ls.integrate(1, 3));
//        System.out.println(ls);
//
//        QuadraticSpline qs2 = newQuadraticSpline(x, y);
//        System.out.println(qs2.evaluateAt(2.5));
//        System.out.println(qs2.differentiate(2.5));
//        System.out.println(qs2.integrate(3));
//        System.out.println(qs2.integrate(1, 3));
//        System.out.println(qs2);
//
//        CubicSpline cs = CubicSpline.newNotAKnotSpline(x, y);
//        System.out.println(cs.evaluateAt(2.5));
//        System.out.println(cs.differentiate(2.5));
//        System.out.println(cs.integrate(3));
//        System.out.println(cs.integrate(1, 3));
//        System.out.println(cs);

        x = new double[]{0.0, 10.0, 15.0, 20.0, 22.5, 30.0};
        y = new double[]{0.0, 227.04, 362.78, 517.35, 602.97, 901.67};

        //
        LinearSpline ls = LinearSpline.newLinearSpline(x, y);
        System.out.println(ls.evaluateAt(16));
        System.out.println(ls.differentiate(16));
        System.out.println(ls.integrate(11, 16));
        System.out.println(ls);

        QuadraticSpline qs2 = newNaturalSpline(x, y);
        System.out.println(qs2.evaluateAt(16));
        System.out.println(qs2.differentiate(16));
        System.out.println(qs2.integrate(11, 16));
        System.out.println(qs2);

        CubicSpline cs = CubicSpline.newNotAKnotSpline(x, y);
        System.out.println(cs.evaluateAt(16));
        System.out.println(cs.differentiate(16));
        System.out.println(cs.integrate(11, 16));
        System.out.println(cs);

        double[] yi = new double[31];
        for(int i = 0; i <= 30; ++i) {
            yi[i] = qs2.evaluateAt(i);
        }

        System.out.println(Arrays.toString(yi));
    }
}
