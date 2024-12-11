package com.wildbitsfoundry.etk4j.control;

import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.linearalgebra.EigenvalueDecompositionDense;
import com.wildbitsfoundry.etk4j.math.linearalgebra.MatrixDense;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;
import com.wildbitsfoundry.etk4j.util.DoubleArrays;

import java.util.Arrays;

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
     * @throws IllegalArgumentException     If the length of the initial conditions is different from the number of states.
     * @see <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lsim.html">lsim</a>
     */
    protected TimeResponse lSim(double[][] input, double[] time, double[] initialConditions,
                                StateSpace ss, IntegrationMethod integrationMethod) {
        double[][] U = DoubleArrays.transpose(input);

        if (U.length != time.length) {
            throw new IllegalArgumentException("The input array and the time array must have the same length.");
        }
        if (time.length == 0) {
            throw new IllegalArgumentException("The time array must have at least one element.");
        }

        MatrixDense A = ss.getA();
        MatrixDense B = ss.getB();
        MatrixDense C = ss.getC();
        MatrixDense D = ss.getD();

        final int noStates = A.getRowCount();
        final int noInputs = B.getColumnCount();
        final int noSteps = time.length;

        // initial conditions
        double[] x0 = initialConditions == null ? new double[noStates] : initialConditions;
        double[][] xOut = new double[noSteps][noStates];

        if (x0.length != noStates) {
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
            if (!DoubleArrays.allClose(delta, 1.0)) {
                throw new NonUniformTimeStepsException("Only uniform time steps are supported.");
            }

            switch (integrationMethod) {
                case ZERO_ORDER_HOLD: {
                    A.multiplyEquals(dt);
                    B.multiplyEquals(dt);
                    double[][] M = new double[noStates + noInputs][];
                    for (int i = 0; i < noStates; ++i) {
                        M[i] = DoubleArrays.concatenate(A.getRow(i), B.getRow(i));
                    }
                    for (int i = noStates; i < noStates + noInputs; ++i) {
                        M[i] = new double[noStates + noInputs];
                    }

                    MatrixDense expMT = new MatrixDense(M).transpose().expm();
                    double[][] Ad = expMT.subMatrix(0, noStates - 1, 0, noStates - 1).getAs2DArray();
                    double[][] Bd = expMT.subMatrix(noStates, expMT.getRowCount() - 1, 0, noStates - 1).getAs2DArray();
                    for (int i = 1; i < noSteps; ++i) {
                        xOut[i] = DoubleArrays.addElementWise(dot(xOut[i - 1], Ad), dot(U[i - 1], Bd));
                    }
                    break;
                }
                case INTERPOLATION: {
                    A.multiplyEquals(dt);
                    B.multiplyEquals(dt);
                    double[][] M = new double[noStates + 2 * noInputs][];
                    for (int i = 0; i < noStates; ++i) {
                        M[i] = DoubleArrays.concatenateAll(A.getRow(i), B.getRow(i), new double[noInputs]);
                    }
                    double[][] identity = MatrixDense.Factory.identity(noInputs).getAs2DArray();
                    for (int i = noStates, j = 0; i < noStates + noInputs; ++i, ++j) {
                        M[i] = DoubleArrays.concatenate(new double[noStates + noInputs], identity[j]);
                    }
                    for (int i = noStates + noInputs; i < noStates + 2 * noInputs; ++i) {
                        M[i] = new double[noStates + 2 * noInputs];
                    }

                    MatrixDense expMT = new MatrixDense(M).transpose().expm();
                    double[][] Ad = expMT.subMatrix(0, noStates - 1, 0, noStates - 1).getAs2DArray();
                    double[][] Bd1 = expMT.subMatrix(noStates + noInputs, expMT.getRowCount() - 1, 0, noStates - 1).getAs2DArray();
                    double[][] Bd0 = expMT.subMatrix(noStates, noStates + noInputs - 1, 0, noStates - 1).getAs2DArray();
                    for (int i = 0; i < Bd0.length; ++i) {
                        DoubleArrays.subtractElementWiseInPlace(Bd0[i], Bd1[i]);
                    }
                    for (int i = 1; i < noSteps; ++i) {
                        xOut[i] = DoubleArrays.addElementWise(dot(xOut[i - 1], Ad), dot(U[i - 1], Bd0));
                        DoubleArrays.addElementWiseInPlace(xOut[i], dot(U[i], Bd1));
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
            DoubleArrays.addElementWiseInPlace(yOut[i], dot(U[i], d));
        }
        return new TimeResponse(time, DoubleArrays.transpose(yOut), xOut);
    }

    private static double[] dot(double[] a, double[][] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("The number of elements in a must match the number of rows in b.");
        }
        MatrixDense A = new MatrixDense(a, 1);
        A.multiplyEquals(new MatrixDense(b));
        return A.getArray();
    }

    /**
     * Step response of the continuous-time system with zero initial conditions and 100 default (calculated) time points.
     * The step method assumes that the underlying system is a SISO (Single-Input, Single-Output) system so for a MIMO
     * (Multiple-Input, Multiple-Output) system, the step is calculated for the first input to first output.
     *
     * @return The step response of the system.
     * @see <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.step.html">step</a>
     */
    public StepResponse step() {
        return step(100);
    }

    /**
     * Step response of the continuous-time system with zero initial conditions and the given {@code numberOfpoints}
     * default (calculated) time points. The step method assumes that the underlying system is a SISO
     * (Single-Input, Single-Output) system so for a MIMO (Multiple-Input, Multiple-Output) system, the step is
     * calculated for the first input to first output.
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
     * default (calculated) time points. The step method assumes that the underlying system is a SISO
     * (Single-Input, Single-Output) system so for a MIMO (Multiple-Input, Multiple-Output) system, the step is
     * calculated for the first input to first output.
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
     * Step response of the continuous-time system with zero initial conditions and the given time points. The step
     * method assumes that the underlying system is a SISO (Single-Input, Single-Output) system so for a MIMO
     * (Multiple-Input, Multiple-Output) system, the step is calculated for the first input to first output.
     *
     * @param time The times at which to evaluate the step response.
     * @return The step response of the system.
     * @see <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.step.html">step</a>
     */
    public StepResponse step(double[] time) {
        return step(time, null, time.length);
    }

    /**
     * Step response of the continuous-time system with the given initial conditions and the given time points. The step
     * method assumes that the underlying system is a SISO (Single-Input, Single-Output) system so for a MIMO
     * (Multiple-Input, Multiple-Output) system, the step is calculated for the first input to first output.
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
     *
     * @param time              The time vector. Can be null then
     *                          {@link #generateDefaultResponseTimes(MatrixDense, int)} will be used
     * @param initialConditions Initial conditions of the system.
     * @param numberOfPoints    Number of points to be used in case the time vector is null.
     * @return The StepReponse of the system.
     */
    protected StepResponse step(double[] time, double[] initialConditions, int numberOfPoints) {
        StateSpace ss = this.toStateSpace();
        time = time == null ? generateDefaultResponseTimes(ss.getA(), numberOfPoints) : time;
        double[][] U = new double[1][];
        U[0] = DoubleArrays.ones(time.length);
        TimeResponse lSim = lSim(U, time, initialConditions, ss, IntegrationMethod.ZERO_ORDER_HOLD);
        return new StepResponse(lSim.getTime(), lSim.getResponse()[0]);
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
     */

    /**
     * Default times for the system response.
     *
     * @param A              The input matrix of the State Space system.
     * @param numberOfPoints The number of points to generate.
     * @return The default response times.
     */
    protected static double[] generateDefaultResponseTimes(MatrixDense A, int numberOfPoints) {
        EigenvalueDecompositionDense eig = A.eig();
        double[] realEig = eig.getRealEigenvalues();
        for (int i = 0; i < realEig.length; ++i) {
            realEig[i] = Math.abs(realEig[i]);
        }
        double r = DoubleArrays.min(realEig);
        if (r == 0.0) {
            r = 1.0;
        }
        double tc = 1.0 / r;
        return DoubleArrays.linSpace(0.0, 7 * tc, numberOfPoints);
    }

    /**
     * Frequency response of the system. The default value for the number of points is 200. For MIMO systems this method
     * calculates the response from the first input to the first output.
     *
     * @return The complex frequency response of the system and the frequency at which each point was calculated.
     */
    public FrequencyResponse calculateFrequencyResponse() {
        return this.calculateFrequencyResponse(200);
    }

    /**
     * Frequency response of the system.For MIMO systems this method calculates the response from the first input to the
     * first output.
     *
     * @param numberOfPoints The number of points in which to evaluate the system.
     * @return The complex frequency response of the system and the frequency at which each point was calculated.
     */
    public FrequencyResponse calculateFrequencyResponse(int numberOfPoints) {
        double[] frequencies = findFrequencies(numberOfPoints);
        return this.calculateFrequencyResponse(frequencies);
    }

    /**
     * Evaluate the system at a given frequency. For MIMO systems, this method calculates the response from the first
     * input to the first output.
     *
     * @param w The frequency at which to evaluate the system.
     * @return The complex frequency response of the system.
     */
    public abstract Complex evaluateAt(double w);

    /**
     * Magnitude of the system. For MIMO systems, this method calculates the magnitude from the first
     * input to the first output.
     * @param w Argument at which to evaluate the function.
     * @return The absolute value of the complex response.
     */
    public double calculateMagnitudeAt(double w) {
        return this.evaluateAt(w).abs();
    }

    /**
     * Magnitude of the system. For MIMO systems, this method calculates the magnitude from the first
     * input to the first output.
     * @param w Argument at which to evaluate the function.
     * @return The absolute value of the complex response.
     */
    public double[] calculateMagnitudeAt(double[] w) {
        double[] magnitude = new double[w.length];
        for (int i = 0; i < w.length; ++i) {
            magnitude[i] = this.calculateMagnitudeAt(w[i]);
        }
        return magnitude;
    }

    /***
     * Calculate the system wrapped phase response. For MIMO systems, this method calculates the phase from the first
     * input to the first output.
     * @param w the frequencies where the phase needs to be calculated at.
     * @return The phase response of the system in rad/s.
     */
    public double calculatePhaseAt(double w) {
        return this.evaluateAt(w).arg();
    }

    /***
     * Calculate the system wrapped phase response. For MIMO systems, this method calculates the phase from the first
     * input to the first output.
     * @param w the frequencies where the phase needs to be calculated at.
     * @return The phase response of the system in rad/s.
     */
    public double[] calculatePhaseAt(double[] w) {
        double[] phase = new double[w.length];
        for (int i = 0; i < w.length; ++i) {
            phase[i] = this.calculatePhaseAt(w[i]);
        }
        return phase;
    }

    /***
     * Calculate the system wrapped phase response. For MIMO systems, this method calculates the magnitude from the first
     * input to the first output.
     * @param w the frequency where the phase needs to be calculated at.
     * @return The phase response of the system in degrees.
     */
    public double calculatePhaseInDegreesAt(double w) {
        return Math.toDegrees(this.evaluateAt(w).arg());
    }

    /***
     * Calculate the system wrapped phase response. For MIMO systems, this method calculates the magnitude from the first
     * input to the first output.
     * @param w the frequencies where the phase needs to be calculated at.
     * @return The phase response of the system in degrees.
     */
    public double[] calculatePhaseInDegreesAt(double[] w) {
        double[] phase = new double[w.length];
        for (int i = 0; i < w.length; ++i) {
            phase[i] = this.calculatePhaseInDegreesAt(w[i]);
        }
        return phase;
    }

    /**
     * Frequency response of the system. For MIMO systems, this method calculates the response from the first input to
     * the first output.
     *
     * @param w The frequencies at which to evaluate the system.
     * @return The frequency response of the system at each given frequency.
     */
    public FrequencyResponse calculateFrequencyResponse(double[] w) {
        Complex[] response = new Complex[w.length];
        for (int i = 0; i < w.length; ++i) {
            response[i] = this.evaluateAt(w[i]);
        }
        return new FrequencyResponse(response, w);
    }

    /**
     * Bode response of the system. For MIMO systems, this method calculates the response from the first input to the
     * first output.
     *
     * @return The magnitude in dB and phase in degrees of the system.
     */
    public BodeResponse calculateBode() {
        return this.calculateBode(200);
    }

    /**
     * Bode response of the system. For MIMO systems, this method calculates the response from the first input to the
     * first output.
     *
     * @param numberOfPoints The number of points in which to evaluate the system.
     * @return The magnitude in dB and phase in degrees of the system.
     */
    public BodeResponse calculateBode(int numberOfPoints) {
        double[] frequencies = findFrequencies(numberOfPoints);
        return this.calculateBode(frequencies);
    }

    /**
     * Bode response of the system. For MIMO systems, this method calculates the response from the first input to the
     * first output.
     *
     * @param w The frequencies at which to evaluate the system.
     * @return The magnitude in dB and phase in degrees of the system.
     */
    public BodeResponse calculateBode(double[] w) {
        double[] magnitudeIndB = new double[w.length];
        double[] phaseInDegrees = new double[w.length];
        for (int i = 0; i < w.length; ++i) {
            Complex response = this.evaluateAt(w[i]);
            magnitudeIndB[i] = 20 * Math.log10(response.abs());
            phaseInDegrees[i] = Math.toDegrees(response.arg());
        }
        unwrapPhase(phaseInDegrees);
        return new BodeResponse(magnitudeIndB, phaseInDegrees, w);
    }

    /**
     * Helper method to find the "interesting frequencies" where the system is exhibiting change.
     *
     * @param numberOfPoints The number of points to calculate the frequency at.
     * @return the frequencies where cool things are happening.
     */
    private double[] findFrequencies(int numberOfPoints) {
        ZeroPoleGain zpk = this.toZeroPoleGain();
        Complex[] ep = zpk.getPoles();
        Complex[] tz = zpk.getZeros();

        if (ep.length == 0) {
            ep = new Complex[]{Complex.fromReal(-1000)};
        }
        Complex[] ez = ComplexArrays.concatenate(
                Arrays.stream(ep).filter(c -> c.imag() >= 0.0).toArray(Complex[]::new),
                Arrays.stream(tz).filter(c -> c.abs() < 1e5 && c.imag() >= 0.0).toArray(Complex[]::new)
        );
        int[] integ = Arrays.stream(ez).mapToInt(c -> c.abs() < 1e-8 ? 1 : 0).toArray();
        double[] argument = new double[ez.length];
        for (int i = 0; i < argument.length; ++i) {
            argument[i] = 3.0 * Math.abs(ez[i].real() + integ[i]) + 1.5 * ez[i].imag();
        }
        int hiFreq = roundFrequency(Math.log10(DoubleArrays.max(argument)) + 0.5);

        for (int i = 0; i < argument.length; ++i) {
            argument[i] = Math.abs(ez[i].add(integ[i]).real()) + 2 * ez[i].imag();
        }
        int loFreq = roundFrequency(Math.log10(0.1 * DoubleArrays.min(argument)) - 0.5);

        return DoubleArrays.logSpace(loFreq, hiFreq, numberOfPoints);
    }

    private static int roundFrequency(double d) {
        // get numbers after the decimal point
        double decimal = d - Math.floor(d);
        if (Math.abs(decimal) == 0.5) {
            return (int) MathETK.roundEven(d);
        } else {
            return (int) Math.round(d);
        }
    }

    /**
     * Single-input Single-output system time response. For MIMO systems, the response is calculated for the first input
     * to first output.
     * @param input The values of the input vs time.
     * @param time The array of time points.
     * @return The time domain response of the system.
     */
    public SISOTimeResponse simulateTimeResponse(double[] input, double[] time) {
        return simulateTimeResponse(input, time, IntegrationMethod.INTERPOLATION);
    }

    /**
     * Single-input Single-output system time response. For MIMO systems, the response is calculated for the first input
     * to first output.
     * @param input The values of the input vs time.
     * @param time The array of time points.
     * @param initialConditions The initial conditions of the system.
     * @return The time domain response of the system.
     */
    public SISOTimeResponse simulateTimeResponse(double[] input, double[] time,
                                                 double[] initialConditions) {
        return simulateTimeResponse(input, time, initialConditions, IntegrationMethod.INTERPOLATION);
    }

    /**
     * Single-input Single-output system time response. For MIMO systems, the response is calculated for the first input
     * to first output.
     * @param input The values of the input vs time.
     * @param time The array of time points.
     * @param integrationMethod The {@link com.wildbitsfoundry.etk4j.control.LinearTimeInvariantSystem.IntegrationMethod}.
     *                          default is INTERPOLATION.
     * @return The time domain response of the system.
     */
    public SISOTimeResponse simulateTimeResponse(double[] input, double[] time,
                                                 IntegrationMethod integrationMethod) {
        return simulateTimeResponse(input, time, null, integrationMethod);
    }

    /**
     * Single-input Single-output system time response. For MIMO systems, the response is calculated for the first input
     * to first output.
     * @param input The values of the input vs time.
     * @param time The array of time points.
     * @param initialConditions The initial conditions of the system.
     * @param integrationMethod The {@link com.wildbitsfoundry.etk4j.control.LinearTimeInvariantSystem.IntegrationMethod}.
     * @return The time domain response of the system.
     */
    public SISOTimeResponse simulateTimeResponse(double[] input, double[] time,
                                                 double[] initialConditions,
                                                 IntegrationMethod integrationMethod) {
        double[][] U = new double[1][time.length];
        U[0] = input;
        TimeResponse tr = lSim(U, time, initialConditions, this.toStateSpace(), integrationMethod);
        return new SISOTimeResponse(tr.getTime(), tr.getResponse()[0], tr.getEvolutionOfStateVector());
    }

    /**
     * Unwraps frequency in degrees. For example, consider a system where the phase response goes from 0°
     * to 360°. If the phase is wrapped, it will only vary between 180° and -180° e.g. (Careful, ASCII art coming your way)
     * <pre>
     *  Expected:                     Actual:
     *      |                           |
     *      |                           |
     *      |                           |
     *      ---------------------       ---------------------
     *    0 |*******                    |******* *******
     *  180 |       *                   |       *
     *  360 |       *******             |
     *  </pre>
     *
     * @param phase The phase to unwrap. This operation is done in plase and the phase array contains the result of
     *              unwrapping the frequency.
     */
    public static void unwrapPhase(double[] phase) {
        int length = phase.length;
        double[] dp = new double[length];
        double[] dps = new double[length];
        double[] C = new double[length];
        double[] cumulativeSum = new double[length];

        double cutoff = 180.0;
        int j;

        // incremental phase variation
        for (j = 0; j < length - 1; j++) {
            dp[j] = phase[j + 1] - phase[j];
        }
        // equivalent phase variation in [-pi, pi]
        for (j = 0; j < length - 1; j++) {
            dps[j] = (dp[j] + 180.0) - Math.floor((dp[j] + 180.0) / (2 * 180.0)) * (2 * 180.0) - 180.0;
        }
        // preserve variation sign for +pi vs. -pi
        for (j = 0; j < length - 1; j++) {
            if ((dps[j] == -180.0) && (dp[j] > 0)) {
                dps[j] = 180.0;
            }
        }
        // incremental phase correction
        for (j = 0; j < length - 1; j++) {
            C[j] = dps[j] - dp[j];
        }
        // Ignore correction when incremental variation is smaller than cutoff
        for (j = 0; j < length - 1; j++) {
            if (Math.abs(dp[j]) < cutoff) {
                C[j] = 0;
            }
        }
        // Find cumulative sum of deltas
        cumulativeSum[0] = C[0];
        for (j = 1; j < length - 1; j++) {
            cumulativeSum[j] = cumulativeSum[j - 1] + C[j];
        }
        // Integrate corrections and add to P to produce smoothed phase values
        for (j = 1; j < length; j++) {
            phase[j] += cumulativeSum[j - 1];
        }
    }
}
