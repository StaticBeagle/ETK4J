package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

public class ComplexCholeskyDecompositionDense extends ComplexCholeskyDecomposition<ComplexMatrixDense> {

    private ComplexMatrixDense R;
    boolean isSPD;

    public ComplexCholeskyDecompositionDense(ComplexMatrixDense matrix) {
        super(matrix);
        Complex[] A = matrix.getArrayCopy();
        Complex[] ARef = matrix.getArray();
        final int m = matrix.getRowCount();
        final int n = matrix.getColumnCount();

        isSPD = (n == m);
        for(int k = 0; k < n; k++) {
            isSPD = isSPD & A[k * n + k].real() > 0 & A[k * n + k].imag() == 0;
            A[k * n + k] = Complex.fromReal(Math.sqrt(A[k * n + k].real()));
            double mu = 1.0 / A[k * n + k].real();
            for(int j = k + 1; j < n; j++) {
                isSPD = isSPD & ARef[k * n + j].equals(ARef[j * n + k].conj());
                A[k * n + j] = A[k * n + j].multiply(mu);
            }
            for(int i = k + 1; i < n; i++) {
                for(int j = i; j < n; j++) {
                    double real = A[i * n + j].real()
                            - A[k * n + i].real() * A[k * n + j].real() - A[k * n + i].imag() * A[k * n + j].imag();
                    double imag = A[i * n + j].imag()
                            - A[k * n + i].real() * A[k * n + j].imag() + A[k * n + i].imag() * A[k * n + j].real();
                    A[i * n + j] = new Complex(real, imag);
                }
                A[i * n + i] = new Complex(A[i * n + i].real(), 0.0);
            }
            for(int j = 0; j < k; j++) {
                A[k * n + j] = new Complex();
            }
        }
        R = new ComplexMatrixDense(A, n, n);
    }

    public ComplexMatrixDense getR() {
        return R.copy();
    }

    public boolean isSPD() {
        return isSPD;
    }

    public static void main(String[] args) {
        ComplexMatrixDense Arg = ComplexMatrixDense.from2DArray(new Complex[][]{
                {new Complex(0.9767407521087846, 0.4711470863081668), new Complex(0.2794308404798361, 0.4587985898568102), new Complex(0.1685348957730604, 0.2220415588608931), new Complex(0.8646965674039528, 0.5496275298999502)},
                {new Complex(0.7198733467511716, 0.9814006625666408), new Complex(0.01506799866339292, 0.666774973911552), new Complex(0.750545968546426, 0.4842583141948102), new Complex(0.8173900466200722, 0.9463910637201288)},
                {new Complex(0.3444764277139777, 0.7365028309538496), new Complex(0.8429046576872399, 0.1345623113997281), new Complex(0.6191817846917982, 0.174500592061623), new Complex(0.1931980127831494, 0.3833402493069025)},
                {new Complex(0.2529620785541479, 0.8308769763271352), new Complex(0.7385598627121508, 0.850513285400206), new Complex(0.3669833799842998, 0.3169573549259639), new Complex(0.4630012679866645, 0.4444797942481419)}
        });

        // Making it symmetric
        System.out.println("Original Matrix");
        Arg.multiplyEquals(Arg.conjugateTranspose());
        System.out.println(Arg);

        ComplexCholeskyDecompositionDense choleskyDecompositionDense = new ComplexCholeskyDecompositionDense(Arg);
        ComplexMatrixDense R = choleskyDecompositionDense.getR();
        System.out.println("R");
        System.out.println(R);
        System.out.println("RH * R");
        System.out.println(R.conjugateTranspose().multiply(R));
        System.out.println("Is symmetric");
        System.out.println(choleskyDecompositionDense.isSPD());
    }
}
