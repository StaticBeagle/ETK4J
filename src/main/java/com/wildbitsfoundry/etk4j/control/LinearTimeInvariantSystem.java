package com.wildbitsfoundry.etk4j.control;

import com.wildbitsfoundry.etk4j.math.linearalgebra.EigenvalueDecomposition;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrices;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;
import com.wildbitsfoundry.etk4j.util.NumArrays;

/**
 * The {@code LinearTimeInvariantSystem} represents and LTI system and provides methods to simulate the time response
 * of the given system.
 */
public abstract class LinearTimeInvariantSystem {

    public enum IntegrationMethod {
        ZERO_ORDER_HOLD,
        INTERPOLATION
    }

    public abstract StateSpace toStateSpace();

    public abstract TransferFunction toTransferFunction();

    public abstract ZeroPoleGain toZeroPoleGain();

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
     */

    /**
     * Simulate time response of a continuous-time system.
     *
     * @param input             Array describing the input at every time step. For multiple inputs, each row of this
     *                          array represents an input to the system.
     * @param time              The time at which to evaluate the system.
     * @param initialConditions Initial conditions of the system.
     * @param ss                State Space representation of the system.
     * @param integrationMethod Integration method.
     * @return The {@link TimeResponse} Of the system.
     * @throws IllegalArgumentException     If the length of the input array does not match the length of the time array.
     * @throws IllegalArgumentException     If the time array is empty.
     * @throws IllegalArgumentException     If the initial time is negative.
     * @throws NonUniformTimeStepsException If the step of the time array ore not uniform.
     * @throws IllegalArgumentException If the length of the initial conditions is different from the number of states.
     * @see <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lsim.html">lsim</a>
     */
    protected TimeResponse lsim(double[][] input, double[] time, double[] initialConditions,
                                StateSpace ss, IntegrationMethod integrationMethod) {
        double[][] U = NumArrays.transpose(input);

        if (U.length != time.length) {
            throw new IllegalArgumentException("The input array and the time array must have the same length.");
        }
        if (time.length == 0) {
            throw new IllegalArgumentException("The time array must have at least one element.");
        }

        Matrix A = ss.getA();
        Matrix B = ss.getB();
        Matrix C = ss.getC();
        Matrix D = ss.getD();

        final int noStates = A.getRowCount();
        final int noInputs = B.getColumnCount();
        final int noSteps = time.length;

        // initial conditions
        double[] x0 = initialConditions == null ? new double[noStates] : initialConditions;
        double[][] xOut = new double[noSteps][noStates];

        if(x0.length != noStates) {
            throw new IllegalArgumentException("The number of initial conditions is different from the number of states.");
        }
        if (time[0] == 0.0) {
            xOut[0] = x0;
        } else if (time[0] > 0.0) {
            xOut[0] = dot(x0, A.transpose().multiply(time[0]).expm().getAs2DArray());
        } else {
            throw new IllegalArgumentException("Initial time must be non negative.");
        }

        if (noSteps > 1) {
            double dt = time[1] - time[0];
            double[] delta = new double[time.length - 2];
            for (int i = 1; i < time.length - 1; ++i) {
                delta[i - 1] = (time[i + 1] - time[i]) / dt;
            }
            if (!NumArrays.allClose(delta, 1.0)) {
                throw new NonUniformTimeStepsException("Only uniform time steps are supported.");
            }

            switch (integrationMethod) {
                case ZERO_ORDER_HOLD: {
                    A.multiplyEquals(dt);
                    B.multiplyEquals(dt);
                    double[][] M = new double[noStates + noInputs][];
                    for (int i = 0; i < noStates; ++i) {
                        M[i] = NumArrays.concatenate(A.getRow(i), B.getRow(i));
                    }
                    for (int i = noStates; i < noStates + noInputs; ++i) {
                        M[i] = new double[noStates + noInputs];
                    }

                    Matrix expMT = new Matrix(M).transpose().expm();
                    double[][] Ad = expMT.subMatrix(0, noStates - 1, 0, noStates - 1).getAs2DArray();
                    double[][] Bd = expMT.subMatrix(noStates, expMT.getRowCount() - 1, 0, noStates - 1).getAs2DArray();
                    for (int i = 1; i < noSteps; ++i) {
                        xOut[i] = NumArrays.add(dot(xOut[i - 1], Ad), dot(U[i - 1], Bd));
                    }
                    break;
                }
                case INTERPOLATION: {
                    A.multiplyEquals(dt);
                    B.multiplyEquals(dt);
                    double[][] M = new double[noStates + 2 * noInputs][];
                    for (int i = 0; i < noStates; ++i) {
                        M[i] = NumArrays.concatenateAll(A.getRow(i), B.getRow(i), new double[noInputs]);
                    }
                    double[][] identity = Matrices.identity(noInputs).getAs2DArray();
                    for (int i = noStates, j = 0; i < noStates + noInputs; ++i, ++j) {
                        M[i] = NumArrays.concatenate(new double[noStates + noInputs], identity[j]);
                    }
                    for (int i = noStates + noInputs; i < noStates + 2 * noInputs; ++i) {
                        M[i] = new double[noStates + 2 * noInputs];
                    }

                    Matrix expMT = new Matrix(M).transpose().expm();
                    double[][] Ad = expMT.subMatrix(0, noStates - 1, 0, noStates - 1).getAs2DArray();
                    double[][] Bd1 = expMT.subMatrix(noStates + noInputs, expMT.getRowCount() - 1, 0, noStates - 1).getAs2DArray();
                    double[][] Bd0 = expMT.subMatrix(noStates, noStates + noInputs - 1, 0, noStates - 1).getAs2DArray();
                    for (int i = 0; i < Bd0.length; ++i) {
                        NumArrays.subtractElementWiseInPlace(Bd0[i], Bd1[i]);
                    }
                    for (int i = 1; i < noSteps; ++i) {
                        xOut[i] = NumArrays.add(dot(xOut[i - 1], Ad), dot(U[i - 1], Bd0));
                        NumArrays.addElementWiseInPlace(xOut[i], dot(U[i], Bd1));
                    }
                    break;
                }
                default:
                    throw new IllegalArgumentException("Unknown integration method.");
            }
        }
        double[][] yOut = new double[noSteps][noStates];
        double[][] c = C.transpose().getAs2DArray();
        double[][] d = D.transpose().getAs2DArray();
        for (int i = 0; i < noSteps; ++i) {
            yOut[i] = dot(xOut[i], c);
            NumArrays.addElementWiseInPlace(yOut[i], dot(U[i], d));
        }
        return new TimeResponse(time, NumArrays.transpose(yOut), xOut);
    }

    private static double[] dot(double[] a, double[][] b) {
        if(a.length != b.length) {
            throw new IllegalArgumentException("The number of elements in a must match the number of rows in b.");
        }
        Matrix A = new Matrix(a, 1);
        A.multiplyEquals(new Matrix(b));
        return A.getArray();
    }

    /**
     * Step response of the continuous-time system with zero initial conditions and 100 default (calculated) time points.
     *
     * @return The step response of the system.
     * @see <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.step.html">step</a>
     */
    public StepResponse step() {
        return step(100);
    }

    /**
     * Step response of the continuous-time system with zero initial conditions and the given {@code numberOfpoints}
     * default (calculated) time points.
     *
     * @param numberOfPoints The number of points in which to evaluate the step response.
     * @return The step response of the system.
     * @see <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.step.html">step</a>
     */
    public StepResponse step(int numberOfPoints) {
        return step(null, null, numberOfPoints);
    }

    /**
     * Step response of the continuous-time system with the given initial conditions and the given {@code numberOfpoints}
     * default (calculated) time points.
     *
     * @param initialConditions The initial conditions of the system.
     * @param numberOfPoints    The number of points in which to evaluate the step response.
     * @return The step response of the system.
     * @see <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.step.html">step</a>
     */
    public StepResponse step(double[] initialConditions, int numberOfPoints) {
        return step(null, initialConditions, numberOfPoints);
    }

    /**
     * Step response of the continuous-time system with zero initial conditions and the given time points.
     *
     * @param time The times at which to evaluate the step response.
     * @return The step response of the system.
     * @see <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.step.html">step</a>
     */
    public StepResponse step(double[] time) {
        return step(time, null, time.length);
    }

    /**
     * Step response of the continuous-time system with the given initial conditions and the given time points.
     *
     * @param time              The times at which to evaluate the step response.
     * @param initialConditions The initial conditions of the system.
     * @return The step response of the system.
     * @see <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.step.html">step</a>
     */
    public StepResponse step(double[] time, double[] initialConditions) {
        return step(time, initialConditions, time.length);
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
     */

    /**
     * Helper function to support all the step overloads.
     * @param time The time vector. Can be null then
     * {@link #generateDefaultResponseTimes(Matrix, int)} will be used
     * @param initialConditions Initial conditions of the system.
     * @param numberOfPoints Number of points to be used in case the time vector is null.
     * @return The StepReponse of the system.
     */
    protected StepResponse step(double[] time, double[] initialConditions, int numberOfPoints) {
        StateSpace ss = this.toStateSpace();
        time = time == null ? generateDefaultResponseTimes(ss.getA(), numberOfPoints) : time;
        double[][] U = new double[1][];
        U[0] = NumArrays.ones(time.length);
        TimeResponse lSim = lsim(U, time, initialConditions, ss, IntegrationMethod.ZERO_ORDER_HOLD);
        return new StepResponse(lSim.getTime(), lSim.getResponse()[0]);
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
     */

    /**
     * Default times for the system response.
     * @param A The input matrix of the State Space system.
     * @param numberOfPoints The number of points to generate.
     * @return The default response times.
     */
    protected static double[] generateDefaultResponseTimes(Matrix A, int numberOfPoints) {
        EigenvalueDecomposition eig = A.eig();
        double[] realEig = eig.getRealEigenvalues();
        for (int i = 0; i < realEig.length; ++i) {
            realEig[i] = Math.abs(realEig[i]);
        }
        double r = NumArrays.min(realEig);
        if (r == 0.0) {
            r = 1.0;
        }
        double tc = 1.0 / r;
        return NumArrays.linSpace(0.0, 7 * tc, numberOfPoints);
    }
}
