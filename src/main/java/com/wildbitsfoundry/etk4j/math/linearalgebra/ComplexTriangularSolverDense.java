package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

import static com.wildbitsfoundry.etk4j.util.ComplexArrays.zeros;

public class ComplexTriangularSolverDense {
    private ComplexTriangularSolverDense() {}

    public static ComplexMatrixDense backSubstitutionSolve(ComplexMatrixDense R, ComplexMatrixDense B) {
        int n = R.getRowCount();
        int m = B.getColumnCount();

        Complex[][] X = zeros(n, m);
        for (int col = 0; col < m; col++) {
            for (int i = n - 1; i >= 0; i--) {
                X[i][col] = B.get(i, col);
                for (int j = i + 1; j < n; j++) {
                    X[i][col].subtractEquals(R.get(i, j).multiply(X[j][col]));
                }
                X[i][col].divideEquals(R.get(i, i));
            }
        }
        return new ComplexMatrixDense(X);
    }

    public static ComplexMatrixDense forwardSubstitutionSolve(ComplexMatrixDense L, ComplexMatrixDense B) {
        int n = L.getRowCount();
        int m = B.getColumnCount();

        Complex[][] X = zeros(n, m);
        for (int col = 0; col < m; col++) {
            for (int i = 0; i < n; i++) {
                X[i][col] = B.get(i, col);
                for (int j = 0; j < i; j++) {
                    X[i][col].subtractEquals(L.get(i, j).multiply(X[j][col]));
                }
                X[i][col].divideEquals(L.get(i, i));
            }
        }
        return new ComplexMatrixDense(X);
    }
}
