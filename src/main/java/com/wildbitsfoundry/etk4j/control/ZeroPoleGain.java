package com.wildbitsfoundry.etk4j.control;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;

import static com.wildbitsfoundry.etk4j.util.ComplexArrays.deepCopy;

/**
 * The {@code ZeroPoleGain} (ZPK) class represents a Linear Time Invariant system broken into its zeros, poles, and gain.
 */
public class ZeroPoleGain extends LinearTimeInvariantSystem {
    private final Complex[] zeros;
    private final Complex[] poles;
    private final double gain;

    public ZeroPoleGain(Complex[] zeros, Complex[] poles, double gain) {
        this.zeros = deepCopy(zeros);
        this.poles = deepCopy(poles);
        this.gain = gain;
    }

    /**
     * Zeros of the system.
     * @return The zeros of the system.
     */
    public Complex[] getZeros() {
        return deepCopy(zeros);
    }

    /**
     * Poles of the system.
     * @return The poles of the system.
     */
    public Complex[] getPoles() {
        return deepCopy(poles);
    }

    /**
     * Gain of the system.
     * @return The scalar gain of the system.
     */
    public double getGain() {
        return gain;
    }

    /**
     * {@link StateSpace} representation of the system.
     * @return The State-Space representation of the system. This method converts the zpk representation into a
     * {@link TransferFunction} and then calls {@link TransferFunction#toStateSpace()}.
     */
    @Override
    public StateSpace toStateSpace() {
        return new TransferFunction(this).toStateSpace();
    }

    /**
     * {@link TransferFunction} representation of the system.
     * @return The transfer function of the system.
     */
    @Override
    public TransferFunction toTransferFunction() {
        return new TransferFunction(this);
    }

    /**
     * {@link ZeroPoleGain} representation of the system.
     * @return Itself.
     */
    @Override
    public ZeroPoleGain toZeroPoleGain() {
        return this;
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
     */

    /**
     * Evaluate the system at a given frequency.
     *
     * @param w The frequency at which to evaluate the system.
     * @return The complex frequency response of the system.
     */
    @Override
    public Complex evaluateAt(double w) {
        Complex num = Polynomial.polyvalFromRoots(zeros, Complex.fromImaginary(w));
        Complex den = Polynomial.polyvalFromRoots(poles, Complex.fromImaginary(w));
        num.divideEquals(den);
        num.multiplyEquals(gain);
        return num;
    }
}
