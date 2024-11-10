package com.wildbitsfoundry.etk4j.control;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.linearalgebra.ComplexMatrixDense;
import com.wildbitsfoundry.etk4j.math.linearalgebra.MatrixDense;
import com.wildbitsfoundry.etk4j.math.linearalgebra.NonSquareMatrixException;
import com.wildbitsfoundry.etk4j.util.DoubleArrays;

import java.util.Arrays;

/**
 * The {@code StateSpace} class represents a Linear Time Invariant System in state-space form.
 */
public class StateSpace extends LinearTimeInvariantSystem {

    private MatrixDense A;
    private MatrixDense B;
    private MatrixDense C;
    private MatrixDense D;

    /**
     * Constructs a null system where all the matrices are empty matrices.
     */
    public StateSpace() {
        this.A = MatrixDense.empty();
        this.B = MatrixDense.empty();
        this.C = MatrixDense.empty();
        this.D = MatrixDense.empty();
    }

    /**
     * Constructs a {@code StateSpace} system.
     * @param A The state matrix in array form.
     * @param B The input to state matrix in array form.
     * @param C The state to output matrix in array form.
     * @param D The feed through matrix in array form.
     */
    public StateSpace(double[][] A, double[][] B, double[][] C, double[][] D) {
        this.A = A == null ? null : new MatrixDense(A);
        this.B = B == null ? null : new MatrixDense(B);
        this.C = C == null ? null : new MatrixDense(C);
        this.D = D == null ? null : new MatrixDense(D);
    }

    /**
     * Constructs a {@code StateSpace} system.
     * @param A The state matrix.
     * @param B The input to state matrix.
     * @param C The state to output matrix.
     * @param D The feed through matrix.
     */
    public StateSpace(MatrixDense A, MatrixDense B, MatrixDense C, MatrixDense D) {
        this.A = A == null ? null : new MatrixDense(A);
        this.B = B == null ? null : new MatrixDense(B);
        this.C = C == null ? null : new MatrixDense(C);
        this.D = D == null ? null : new MatrixDense(D);
    }

    /**
     * Get the state {@link MatrixDense}.
     * @return A copy of the state matrix A.
     */
    public MatrixDense getA() {
        return new MatrixDense(A);
    }

    /**
     * Get input to state {@link MatrixDense}.
     * @return A copy of the input to state matrix B.
     */
    public MatrixDense getB() {
        return new MatrixDense(B);
    }

    /**
     * Get the state to output {@link MatrixDense}.
     * @return A copy of the state to output matrix C.
     */
    public MatrixDense getC() {
        return new MatrixDense(C);
    }

    /**
     * Get the feed through {@link MatrixDense}.
     * @return A copy of the feed through matrix D.
     */
    public MatrixDense getD() {
        return new MatrixDense(D);
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
        StateSpace ss = normalize(this);
        MatrixDense A = ss.A;
        MatrixDense B = ss.B;
        MatrixDense C = ss.C;
        MatrixDense D = ss.D;
        int nin = D.getColumnCount();
        int nout = D.getRowCount();
        if(input >= nin) {
            throw new IllegalArgumentException("System does not have the input specified.");
        }
        // make SIMO from possibly MIMO system.
        B = B.subMatrix(0, B.getRowCount() - 1, input, input);
        D = D.subMatrix(0, D.getRowCount() - 1, input, input);

        double[][] num;
        double[] den;

        if(!A.isSquared()) {
            throw new NonSquareMatrixException("Matrix A must be a square Matrix.");
        }
        den = A.poly();

        if(B.getRowCount() * B.getColumnCount() == 0 && C.getRowCount() * C.getColumnCount() == 0) {
            if(D.getRowCount() * D.getColumnCount() == 0 && A.getRowCount() * A.getColumnCount() == 0) {
                den = new double[0];
            }
            return new TransferFunction[] { new TransferFunction(D.getArray(), den)};
        }

        int numStates = A.getRowCount();
        num = new double[nout][numStates + 1];
        TransferFunction[] tfs = new TransferFunction[nout];
        for(int k = 0; k < nout; ++k) {
            double[] Ck = C.getRow(k);
            double Dk = D.get(k, 0);
            num[k] = A.subtract(new MatrixDense(dot(B.getArray(), Ck))).poly();
            DoubleArrays.addElementWiseInPlace(num[k], DoubleArrays.multiplyElementWise(den, Dk - 1));
            tfs[k] = new TransferFunction(num[k], den);
        }
        return tfs;
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
     */
    private StateSpace normalize(StateSpace ss) {
        StateSpace result = new StateSpace(ss.A, ss.B, ss.C, ss.D);
        Integer[] tmp = shapeOrNull(result.A);
        Integer Ma = tmp[0];
        Integer Na = tmp[1];

        tmp = shapeOrNull(result.B);
        Integer Mb = tmp[0];
        Integer Nb = tmp[1];

        tmp = shapeOrNull(result.C);
        Integer Mc = tmp[0];
        Integer Nc = tmp[1];

        tmp = shapeOrNull(result.D);
        Integer Md = tmp[0];
        Integer Nd = tmp[1];

        Integer p = choiceNotNull(Ma, Mb, Nc);
        Integer q = choiceNotNull(Nb, Nd);
        Integer r = choiceNotNull(Mc, Md);
        if(p == null || q == null || r == null) {
            throw new RuntimeException("Not enough information on the system.");
        }

        if(result.A == null) {
            result.A = MatrixDense.empty();
        }
        if(result.B == null) {
            result.B = MatrixDense.empty();
        }
        if(result.C == null) {
            result.C = MatrixDense.empty();
        }
        if(result.D == null) {
            result.D = MatrixDense.empty();
        }
        result.A = restore(result.A, p, p);
        result.B = restore(result.B, p, q);
        result.C = restore(result.C, r, p);
        result.D = restore(result.D, r, q);

        return result;
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
     */
    private Integer[] shapeOrNull(MatrixDense m) {
        return m != null ? new Integer[] {m.getRowCount(), m.getColumnCount()} : new Integer[2];
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
     */
    private Integer choiceNotNull(Integer ...params) {
        for(Integer i : params) {
            if(i != null) {
                return i;
            }
        }
        return null;
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
     */
    private MatrixDense restore(MatrixDense m, int rows, int cols) {
        if(m.isEmpty()) {
            return new MatrixDense(rows, cols);
        } else {
            if(m.getRowCount() != rows && m.getColumnCount() != cols) {
                throw new RuntimeException("The input arrays have incompatible shapes.");
            }
            return m;
        }
    }

    private static double[][] dot(double[] a, double[] b) {
        double[][] result = new double[a.length][b.length];
        for(int i = 0; i < a.length; ++i) {
            result[i] = DoubleArrays.multiplyElementWise(b, a[i]);
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
        return lSim(input, time, initialConditions, this, integrationMethod);
    }

    @Override
    public Complex evaluateAt(double w) {
        return evaluateMIMOAt(w)[0];
    }

    /**
     * Evaluate the system at a given frequency.
     * @param w The frequency at which to evaluate the system.
     * @return The complex frequency response of the system.
     */
    public Complex[] evaluateMIMOAt(double w) {
        ComplexMatrixDense inner = MatrixDense.identity(A.getRowCount()).multiply(Complex.fromImaginary(w)).subtract(A).inv();
        ComplexMatrixDense outer = C.multiply(inner).multiply(B);
        return ComplexMatrixDense.fromRealMatrix(D).add(outer).transpose().getArray();
    }

    /**
     * Magnitude of the system.
     * @param w Argument at which to evaluate the function.
     * @return The absolute value of the complex response.
     */
    public double[] calculateMagnitudeMIMOAt(double w) {
        return Arrays.stream(evaluateMIMOAt(w)).mapToDouble(Complex::abs).toArray();
    }

    /**
     * Magnitude of the system.
     * @param w Argument at which to evaluate the function.
     * @return The absolute value of the complex response. Each row in the output represents one I/O combination.
     */
    public double[][] calculateMagnitudeMIMOAt(double[] w) {
        int dim = D.getRowCount() * D.getColumnCount();
        double[][] result = new double[dim][w.length];
        for (int i = 0; i < w.length; ++i) {
            double[] magnitude = this.calculateMagnitudeMIMOAt(w[i]);
            for(int j = 0; j < dim; ++j) {
                result[j][i] = magnitude[j];
            }
        }
        return result;
    }

    /***
     * Calculate the system wrapped phase response.
     * @param w the frequencies where the phase needs to be calculated at.
     * @return The phase response of the system in rad/s.
     */
    public double[] calculatePhaseMIMOAt(double w) {
        return Arrays.stream(evaluateMIMOAt(w)).mapToDouble(Complex::arg).toArray();
    }

    /***
     * Calculate the system wrapped phase response. 
     * @param w the frequencies where the phase needs to be calculated at.
     * @return The phase response of the system in rad/s. Each row in the output represents one I/O combination.
     */
    public double[][] calculatePhaseMIMOAt(double[] w) {
        int dim = D.getRowCount() * D.getColumnCount();
        double[][] result = new double[dim][w.length];
        for (int i = 0; i < w.length; ++i) {
            double[] phase = this.calculatePhaseMIMOAt(w[i]);
            for(int j = 0; j < dim; ++j) {
                result[j][i] = phase[j];
            }
        }
        return result;
    }

    /***
     * Calculate the system wrapped phase response.
     * @param w the frequency where the phase needs to be calculated at.
     * @return The phase response of the system in degrees.
     */
    public double[] calculatePhaseInDegreesMIMOAt(double w) {
        return Arrays.stream(evaluateMIMOAt(w)).mapToDouble(c -> Math.toDegrees(c.arg())).toArray();
    }

    /***
     * Calculate the system wrapped phase response.
     * @param w the frequencies where the phase needs to be calculated at.
     * @return The phase response of the system in degrees. Each row in the output represents one I/O combination.
     */
    public double[][] calculatePhaseInDegreesMIMOAt(double[] w) {
        int dim = D.getRowCount() * D.getColumnCount();
        double[][] result = new double[dim][w.length];
        for (int i = 0; i < w.length; ++i) {
            double[] phase = this.calculatePhaseInDegreesMIMOAt(w[i]);
            for(int j = 0; j < dim; ++j) {
                result[j][i] = phase[j];
            }
        }
        return result;
    }
}
