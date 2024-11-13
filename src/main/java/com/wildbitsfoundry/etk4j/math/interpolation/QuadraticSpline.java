package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.math.linearalgebra.GaussianElimination;
import com.wildbitsfoundry.etk4j.math.linearalgebra.LUDecompositionDense;
import com.wildbitsfoundry.etk4j.math.linearalgebra.MatrixDense;
import com.wildbitsfoundry.etk4j.util.DoubleArrays;

import java.util.Arrays;

import static com.wildbitsfoundry.etk4j.math.linearalgebra.Matrices.solveLDLtTridiagonalSystem;
import static com.wildbitsfoundry.etk4j.math.linearalgebra.Matrices.solveLDUTridiagonalSystem;

public class QuadraticSpline extends Spline {

    private static final double P5 = 0.5, P33 = 1.0 / 3.0;

    protected QuadraticSpline(double[] x, double[] coefficients) {
        super(x, 3);
        this.coefficients = coefficients;
    }

    public static QuadraticSpline newNaturalSpline(double[] x, double[] y) {
        return newNaturalSplineInPlace(Arrays.copyOf(x, x.length), y);
    }

    public static QuadraticSpline newNaturalSplineInPlace(double[] x, double[] y) {
        final int n = x.length;

        MatrixDense A = new MatrixDense(2 * n - 2, 2 * n - 2);
        double[] b = new double[2 * n - 2];

        for (int j = 0; j <= n - 2; j++) {
            double hx = x[j + 1] - x[j];
            A.set(j, j, hx);
            A.set(j, j + n - 1, hx * hx);
            b[j] = y[j + 1] - y[j];
        }

        for (int j = 0; j < n - 2; j++) {
            double hx = x[j + 1] - x[j];
            A.set(j + n - 1, j, 1);
            A.set(j + n - 1, j + 1, -1);
            A.set(j + n - 1, j + n - 1, 2 * hx);
        }

        A.set(2 * n - 4, n - 1, 1); // c[0] = 0
        A.set(2 * n - 3, 2 * n - 3, 1); // c[n - 1] = 0

        double[] solution = GaussianElimination.solve(A, b);
        // compute coefficients
        double[] coefficients = new double[(n - 1) * 3];
        for (int i = 0, j = 0; i < n - 1; ++i, ++j) {
            coefficients[j] = solution[i + n - 1];
            coefficients[++j] = solution[i];
            coefficients[++j] = y[i];
        }
        return new QuadraticSpline(x, coefficients);
    }

    public static QuadraticSpline newClampedSpline(double[] x, double[] y, double m0, double mn) {
        return newClampedSplineInPlace(Arrays.copyOf(x, x.length), y, m0, mn);
    }

    public static QuadraticSpline newClampedSplineInPlace(double[] x, double[] y, double m0, double mn) {
        final int n = x.length;

        MatrixDense A = new MatrixDense(2 * n - 2, 2 * n - 2);
        double[] b = new double[2 * n - 2];

        for (int j = 0; j <= n - 2; j++) {
            double hx = x[j + 1] - x[j];
            A.set(j, j, hx);
            A.set(j, j + n - 1, hx * hx);
            b[j] = y[j + 1] - y[j];
        }

        for (int j = 0; j < n - 2; j++) {
            double hx = x[j + 1] - x[j];
            A.set(j + n - 1, j, 1);
            A.set(j + n - 1, j + 1, -1);
            A.set(j + n - 1, j + n - 1, 2 * hx);
        }

        A.set(2 * n - 4, 0, 1);
        b[2 * n - 4] = m0; // derivative at the start

        A.set(2 * n - 3, n - 2, 1);
        A.set(2 * n - 3, 2 * n - 3, 2 * (x[n - 1] - x[n - 2]));
        b[2 * n - 3] = mn; // derivative at the end

        double[] solution = GaussianElimination.solve(A, b);
        // compute coefficients
        double[] coefficients = new double[(n - 1) * 3];
        for (int i = 0, j = 0; i < n - 1; ++i, ++j) {
            coefficients[j] = solution[i + n - 1];
            coefficients[++j] = solution[i];
            coefficients[++j] = y[i];
        }
        return new QuadraticSpline(x, coefficients);
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
        double[] x = {0.0, 10.0, 15.0, 20.0, 22.5, 30.0};
        double[] y = {0.0, 227.04, 362.78, 517.35, 602.97, 901.67};
//        double[] x = {0, 1, 2, 3};
//        double[] y = {1, Math.exp(1), Math.exp(2), Math.exp(3)};
//        double[] x = {0, 1, 2, 3, 4};
//        double[] y = {1, 2, 0, 2, 1};

//        QuadraticSpline qs2 = QuadraticSpline.newNaturalSpline(x, y);
//
//        System.out.println(qs2.evaluateAt(3));
//
//        CubicSpline cs = CubicSpline.newNaturalSpline(x, y);
//
//        System.out.println(cs.evaluateAt(3));

        QuadraticSpline qs3 = QuadraticSpline.newClampedSpline(x, y, 1, Math.exp(3));

        System.out.println(qs3.evaluateAt(3));

        System.out.println(qs3);


//        CubicSpline cs3 = CubicSpline.newClampedSpline(x, y, 1, Math.exp(3));
//
//        System.out.println(cs3.evaluateAt(3));

//        QuadraticSpline qs3 = QuadraticSpline.newNaturalSpline(x, y);
//
//        System.out.println(qs3.evaluateAt(3));
//
//        System.out.println(qs3);
//
//        QuadraticSpline qs4 = QuadraticSpline.newMeowMix2(x, y);
//
//        System.out.println(qs4.evaluateAt(3));
//
//        System.out.println(qs4);

    }
}
