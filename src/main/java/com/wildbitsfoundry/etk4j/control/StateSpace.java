package com.wildbitsfoundry.etk4j.control;

import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrices;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;

public class StateSpace extends LinearTimeInvariantSystem {

    private Matrix A;
    private Matrix B;
    private Matrix C;
    private Matrix D;

    public StateSpace() {
        this.A = Matrices.empty();
        this.B = Matrices.empty();
        this.C = Matrices.empty();
        this.D = Matrices.empty();
    }

    public StateSpace(double[][] A, double[][] B, double[][] C, double[][] D) {
        this.A = new Matrix(A);
        this.B = new Matrix(B);
        this.C = new Matrix(C);
        this.D = new Matrix(D);
    }

    public StateSpace(Matrix A, Matrix B, Matrix C, Matrix D) {
        this.A = new Matrix(A);
        this.B = new Matrix(B);
        this.C = new Matrix(C);
        this.D = new Matrix(D);
    }

    public Matrix getA() {
        return new Matrix(A);
    }

    public Matrix getB() {
        return new Matrix(B);
    }

    public Matrix getC() {
        return new Matrix(C);
    }

    public Matrix getD() {
        return new Matrix(D);
    }

    @Override
    public String toString() {
        return String.format("A:%n%s%nB:%n%s%nC:%n%s%nD%n%s%n", A, B, C, D);
    }

    @Override
    protected StateSpace toStateSpace() {
        return this;
    }

    // TODO implement this and make them public
    @Override
    protected TransferFunction toTransferFunction() {
        throw null;
    }

    @Override
    protected ZeroPoleGain toZeroPoleGain() {
        throw null;
    }

    public TimeResponse simulateTimeResponse(double[][] input, double[] time) {
        return simulateTimeResponse(input, time, null);
    }

    public TimeResponse simulateTimeResponse(double[][] input, double[] time, double[] initialConditions) {
        return lsim(input, time, initialConditions, this);
    }
}
