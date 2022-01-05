package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;

import java.util.Arrays;

public class QuadraticSpline extends Spline {

    private static final double P5 = 0.5, P33 = 1.0 / 3.0;


    protected QuadraticSpline(double[] x, double[] y) {
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
        _coefs = X.solve(b).getArray();
    }

    public static QuadraticSpline newQuadraticSpline(double[] x, double[] y) {
        return new QuadraticSpline(Arrays.copyOf(x, x.length), y);
    }

    public static QuadraticSpline newQuadraticSplineInPlace(double[] x, double[] y) {
        return new QuadraticSpline(x, y);
    }

    @Override
    public double evaluateSegmentAt(int i, double x) {
        i *= 3;
        return _coefs[i + 2] + x * (_coefs[i + 1] + x * _coefs[i]);
    }

    @Override
    public double evaluateDerivativeAt(int i, double t) {
		double x = t + _x[i];
		i *= 3;
        return _coefs[i + 1] + x * 2.0 * _coefs[i];
    }

    @Override
    public double evaluateAntiDerivativeAt(int i, double t) {
		double x = t + _x[i];
		i *= 3;
        return x * (_coefs[i + 2] + x * (_coefs[i + 1] * P5 + x * _coefs[i] * P33));
    }

    public static void main(String[] args) {
        QuadraticSpline qs2 = newQuadraticSpline(new double[]{0.0, 10.0, 15.0, 20.0, 22.5, 30.0},
                new double[]{0.0, 227.04, 362.78, 517.35, 602.97, 901.67});
        System.out.println(qs2.evaluateAt(16));
		System.out.println(qs2.differentiate(16));
		System.out.println(qs2.integrate(11, 16));
        System.out.println(qs2);
    }
}
