package com.wildbitsfoundry.etk4j.math.polynomials;

import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import com.wildbitsfoundry.etk4j.math.linearalgebra.LeastSquaresSolver;
import com.wildbitsfoundry.etk4j.math.linearalgebra.MatrixDense;

import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkXYDimensions;

public class ChebyshevPolynomial implements UnivariateFunction {

    private final double[] coefficients;

    public ChebyshevPolynomial(double[] x, double[] y, int n) {
        checkXYDimensions(x, y);
        int dim = x.length;
        // Building the coefficient matrix
        MatrixDense A = new MatrixDense(dim, n + 1);
        for(int i = 0; i < dim; i++) {
            for(int j = 0; j <= n; j++) {
                A.unsafeSet(i, j, evaluateChebyshevFirstKind(j, x[i]));
            }
        }
        MatrixDense c = n < y.length ? LeastSquaresSolver.solve(A, y) : A.solve(y);

        double[] coefficients = new double[n + 1];
        for (int i = 0; i <= n; i++) {
            coefficients[i] = c.get(n - i, 0);
        }
        this.coefficients = coefficients;
    }

    @Override
    public double evaluateAt(double x) {
        int n = coefficients.length;
        double c0;
        double c1;
        if(n == 1) {
            c0 = coefficients[0];
            c1 = 0;
        } else if(n == 2) {
            c0 = coefficients[1];
            c1 = coefficients[0];
        } else {
            c0 = coefficients[1];
            c1 = coefficients[0];
            for(int i = 2; i < n; i++) {
                double tmp = c0;
                c0 = coefficients[i] - c1;
                c1 = tmp + 2 * x * c1;
            }
        }
        return c0 + c1 * x;
    }

    public double[] getCoefficients() {
        double[] coefficients = new double[this.coefficients.length];
        System.arraycopy(this.coefficients, 0, coefficients, 0, this.coefficients.length);
        return coefficients;
    }

    public static double[] computeNodes(int n) {
        return computeNodes(-1, 1, n);
    }

    // Compute Chebyshev nodes for a given interval [a, b]
    public static double[] computeNodes(double a, double b, int n) {
        double[] nodes = new double[n];
        for (int k = 0; k < n; k++) {
            nodes[k] = 0.5 * (a + b) + 0.5 * (b - a) * Math.cos((2.0 * k + 1) / (2.0 * n) * Math.PI);
        }
        return nodes;
    }

    public static ChebyshevPolynomial chebyFit(double[] x, double[] y, int n) {
        return new ChebyshevPolynomial(x, y, n);
    }

    // Evaluate the Chebyshev polynomial of the first kind Tn(x)
    public static double evaluateChebyshevFirstKind(int n, double x) {
        if (n == 0) {
            return 1;
        }
        if (n == 1) {
            return x;
        }
        double T0 = 1, T1 = x, Tn = 0;
        for (int i = 2; i <= n; i++) {
            Tn = 2 * x * T1 - T0;
            T0 = T1;
            T1 = Tn;
        }
        return Tn;
    }
}
