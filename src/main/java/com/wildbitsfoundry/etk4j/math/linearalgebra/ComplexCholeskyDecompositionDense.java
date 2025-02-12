package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

import static com.wildbitsfoundry.etk4j.math.linearalgebra.ComplexTriangularSolverDense.backSubstitutionSolve;
import static com.wildbitsfoundry.etk4j.math.linearalgebra.ComplexTriangularSolverDense.forwardSubstitutionSolve;

public class ComplexCholeskyDecompositionDense extends ComplexCholeskyDecomposition<ComplexMatrixDense> {

    private final ComplexMatrixDense R;
    boolean isSPD;

    public ComplexCholeskyDecompositionDense(ComplexMatrixDense matrix) {
        super(matrix);
        Complex[] A = matrix.getArrayCopy();
        Complex[] ARef = matrix.getArray();
        final int m = matrix.getRowCount();
        final int n = matrix.getColumnCount();

        if(m != n) {
            throw new NonSquareMatrixException("Matrix must be squared");
        }

        isSPD = true;
        for (int k = 0; k < n; k++) {
            isSPD = isSPD & A[k * n + k].real() > 0 & A[k * n + k].imag() == 0;
            A[k * n + k] = Complex.fromReal(Math.sqrt(A[k * n + k].real()));
            double mu = 1.0 / A[k * n + k].real();
            for (int j = k + 1; j < n; j++) {
                isSPD = isSPD & ARef[k * n + j].equals(ARef[j * n + k].conj());
                A[k * n + j] = A[k * n + j].multiply(mu);
            }
            for (int i = k + 1; i < n; i++) {
                for (int j = i; j < n; j++) {
                    double real = A[i * n + j].real()
                            - A[k * n + i].real() * A[k * n + j].real() - A[k * n + i].imag() * A[k * n + j].imag();
                    double imag = A[i * n + j].imag()
                            - A[k * n + i].real() * A[k * n + j].imag() + A[k * n + i].imag() * A[k * n + j].real();
                    A[i * n + j] = new Complex(real, imag);
                }
                A[i * n + i] = new Complex(A[i * n + i].real(), 0.0);
            }
            for (int j = 0; j < k; j++) {
                A[k * n + j] = new Complex();
            }
        }
        R = new ComplexMatrixDense(A, n, n);
    }

    public ComplexMatrixDense getR() {
        return R.copy();
    }

    public ComplexMatrixDense getL() {
        return R.conjugateTranspose();
    }

    public boolean isSPD() {
        return isSPD;
    }

    /**
     * Solve A*X = B
     *
     * @param B A Matrix with as many rows as A and any number of columns.
     * @return X so that L*L'*X = B
     * @throws IllegalArgumentException Matrix row dimensions must agree.
     * @throws RuntimeException         Matrix is not symmetric positive definite.
     */

    public ComplexMatrixDense solve(ComplexMatrixDense B) {
        if (B.getRowCount() != rows) {
            throw new IllegalArgumentException("Matrix row dimensions must agree.");
        }
        if (!isSPD) {
            throw new RuntimeException("Matrix is not symmetric positive definite.");
        }

        // Solve RH * Y = B;
        ComplexMatrixDense Y = forwardSubstitutionSolve(R.conjugateTranspose(), B.copy());

        // Solve R * X = Y;
        return backSubstitutionSolve(R, Y);
    }
}
