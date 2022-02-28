package com.wildbitsfoundry.etk4j.control;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.linearalgebra.ComplexMatrix;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrices;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;
import com.wildbitsfoundry.etk4j.math.linearalgebra.NonSquareMatrixException;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;
import com.wildbitsfoundry.etk4j.util.NumArrays;

import java.util.Arrays;

/**
 * The {@code StateSpace} class represents a Linear Time Invariant System in state-space form.
 */
public class StateSpace extends LinearTimeInvariantSystem {

    private Matrix A;
    private Matrix B;
    private Matrix C;
    private Matrix D;

    /**
     * Constructs a null system where all the matrices are empty matrices.
     */
    public StateSpace() {
        this.A = Matrix.empty();
        this.B = Matrix.empty();
        this.C = Matrix.empty();
        this.D = Matrix.empty();
    }

    /**
     * Constructs a {@code StateSpace} system.
     * @param A The state matrix in array form.
     * @param B The input to state matrix in array form.
     * @param C The state to output matrix in array form.
     * @param D The feed through matrix in array form.
     */
    public StateSpace(double[][] A, double[][] B, double[][] C, double[][] D) {
        this.A = new Matrix(A);
        this.B = new Matrix(B);
        this.C = new Matrix(C);
        this.D = new Matrix(D);
    }

    /**
     * Constructs a {@code StateSpace} system.
     * @param A The state matrix.
     * @param B The input to state matrix.
     * @param C The state to output matrix.
     * @param D The feed through matrix.
     */
    public StateSpace(Matrix A, Matrix B, Matrix C, Matrix D) {
        this.A = new Matrix(A);
        this.B = new Matrix(B);
        this.C = new Matrix(C);
        this.D = new Matrix(D);
    }

    /**
     * Get the state {@link Matrix}.
     * @return A copy of the state matrix A.
     */
    public Matrix getA() {
        return new Matrix(A);
    }

    /**
     * Get input to state {@link Matrix}.
     * @return A copy of the input to state matrix B.
     */
    public Matrix getB() {
        return new Matrix(B);
    }

    /**
     * Get the state to output {@link Matrix}.
     * @return A copy of the state to output matrix C.
     */
    public Matrix getC() {
        return new Matrix(C);
    }

    /**
     * Get the feed through {@link Matrix}.
     * @return A copy of the feed through matrix D.
     */
    public Matrix getD() {
        return new Matrix(D);
    }

    @Override
    public String toString() {
        return String.format("A:%n%s%nB:%n%s%nC:%n%s%nD%n%s%n", A, B, C, D);
    }

    @Override
    public StateSpace toStateSpace() {
        return this;
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
     */

    /**
     * {@link TransferFunction} representation of the {@code StateSpace} system.
     * @param input For Multiple Input systems this represents the index of the input to use.
     * @return An array of Transfer Functions, one transfer function per each output of the State Space system.
     */
    public TransferFunction[] toTransferFunction(int input) {
        //TODO normalizeABCD
        int nin = D.getColumnCount();
        int nout = D.getRowCount();
        if(input >= nin) {
            throw new IllegalArgumentException("System does not have the input specified.");
        }
        // make SIMO from possibly MIMO system.
        Matrix Bb = B.subMatrix(0, B.getRowCount() - 1, input, input);
        Matrix Dd = D.subMatrix(0, D.getRowCount() - 1, input, input);

        double[][] num;
        double[] den;

        if(!A.isSquared()) {
            throw new NonSquareMatrixException("Matrix A must be a square Matrix.");
        }
        den = A.poly();

        // TODO test this
        if(Bb.getRowCount() * Bb.getColumnCount() == 0 && C.getRowCount() * C.getColumnCount() == 0) {
            if(Dd.getRowCount() * Dd.getColumnCount() == 0 && A.getRowCount() * A.getColumnCount() == 0) {
                den = new double[0];
            }
            return new TransferFunction[] { new TransferFunction(D.getArray(), den)};
        }

        int numStates = A.getRowCount();
        num = new double[nout][numStates + 1];
        TransferFunction[] tfs = new TransferFunction[nout];
        for(int k = 0; k < nout; ++k) {
            double[] Ck = C.getRow(k);
            double Dk = Dd.get(k, 0);
            num[k] = A.subtract(new Matrix(dot(Bb.getArray(), Ck))).poly();
            NumArrays.addElementWiseInPlace(num[k], NumArrays.multiplyElementWise(den, Dk - 1));
            tfs[k] = new TransferFunction(num[k], den);
        }
        return tfs;
    }

    private static double[][] dot(double[] a, double[] b) {
        double[][] result = new double[a.length][b.length];
        for(int i = 0; i < a.length; ++i) {
            result[i] = NumArrays.multiplyElementWise(b, a[i]);
        }
        return result;
    }

    /**
     * {@link TransferFunction} representation of the {@code StateSpace} system. This method assumes that the
     * system is a Single Input Single Output (SISO) system thus this is similar to calling
     * {@link #toTransferFunction(int)} with input argument equal to zero and returning the first element
     * of the array of returned {@link TransferFunction} array.
     * @return The SISO Transfer Function of the system.
     */
    @Override
    public TransferFunction toTransferFunction() {
        return this.toTransferFunction(0)[0];
    }

    /**
     * {@link ZeroPoleGain} representation of the {@code StateSpace} system.
     * @param input For Multiple Input systems this represents the index of the input to use.
     * @return An array of {@link ZeroPoleGain}, one {@link ZeroPoleGain} per each output of the State Space system.
     */
    public ZeroPoleGain[] toZeroPoleGain(int input) {
        TransferFunction[] tfs = this.toTransferFunction(input);
        return Arrays.stream(tfs).map(TransferFunction::toZeroPoleGain).toArray(ZeroPoleGain[]::new);
    }

    /**
     * {@link ZeroPoleGain} representation of the {@code StateSpace} system. This method assumes that the
     * system is a Single Input Single Output (SISO) system thus this is similar to calling
     * {@link #toZeroPoleGain(int)} with input argument equal to zero and returning the first element
     * of the array of returned {@link ZeroPoleGain} array.
     * @return The SISO Zero, Pole, gain representation of the system.
     */
    @Override
    public ZeroPoleGain toZeroPoleGain() {
        return this.toTransferFunction().toZeroPoleGain();
    }

    /**
     * Simulate time response of a continuous time system with zero initial conditions and interpolation between time steps.
     * @param input Array describing the input at every time step. For multiple inputs, each row of this
     *              array represents an input to the system.
     * @param time The time vector at which to evaluate the response.
     * @return The time response of the system.
     */
    public TimeResponse simulateTimeResponse(double[][] input, double[] time) {
        return simulateTimeResponse(input, time, IntegrationMethod.INTERPOLATION);
    }

    /**
     * Simulate time response of a continuous time system with interpolation between time steps.
     * @param input Array describing the input at every time step. For multiple inputs, each row of this
     *              array represents an input to the system.
     * @param time The time vector at which to evaluate the response.
     * @param initialConditions The initial conditions of the system.
     * @return The time response of the system.
     */
    public TimeResponse simulateTimeResponse(double[][] input, double[] time, double[] initialConditions) {
        return simulateTimeResponse(input, time, initialConditions, IntegrationMethod.INTERPOLATION);
    }

    /**
     * Simulate time response of a continuous time system with zero initial conditions.
     * @param input Array describing the input at every time step. For multiple inputs, each row of this
     *              array represents an input to the system.
     * @param time The time vector at which to evaluate the response.
     * @param integrationMethod The integration method between time points.
     * @return The time response of the system.
     */
    public TimeResponse simulateTimeResponse(double[][] input, double[] time, IntegrationMethod integrationMethod) {
        return simulateTimeResponse(input, time, null, IntegrationMethod.INTERPOLATION);
    }

    /**
     * Simulate time response of a continuous time system.
     * @param input Array describing the input at every time step. For multiple inputs, each row of this
     *              array represents an input to the system.
     * @param time The time vector at which to evaluate the response.
     * @param initialConditions The initial conditions of the system.
     * @param integrationMethod The integration method between time points.
     * @return The time response of the system.
     */
    public TimeResponse simulateTimeResponse(double[][] input, double[] time, double[] initialConditions,
                                             IntegrationMethod integrationMethod) {
        return lsim(input, time, initialConditions, this, integrationMethod);
    }

    /**
     * Evaluate the system at a given frequency.
     * @param w The frequency at which to evaluate the system.
     * @return The complex frequency response of the system.
     */
    public Complex[] evaluateAt(double w) {
        ComplexMatrix inner = Matrix.identity(A.getRowCount()).multiply(Complex.fromImaginary(w)).subtract(A).inv();
        ComplexMatrix outer = C.multiply(inner).multiply(B);
        return ComplexMatrix.fromRealMatrix(D).add(outer).transpose().getArray();
    }
}
