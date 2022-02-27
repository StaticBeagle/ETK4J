package com.wildbitsfoundry.etk4j.control;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;

import static com.wildbitsfoundry.etk4j.util.ComplexArrays.deepCopy;

// TODO document
public class ZeroPoleGain extends LinearTimeInvariantSystem {
    private final Complex[] zeros;
    private final Complex[] poles;
    private final double gain;

    public ZeroPoleGain(Complex[] zeros, Complex[] poles, double gain) {
        this.zeros = deepCopy(zeros);
        this.poles = deepCopy(poles);
        this.gain = gain;
    }

    public Complex[] getZeros() {
        return deepCopy(zeros);
    }

    public Complex[] getPoles() {
        return deepCopy(poles);
    }

    public double getGain() {
        return gain;
    }

    @Override
    public StateSpace toStateSpace() {
        return new TransferFunction(this).toStateSpace();
    }

    @Override
    public TransferFunction toTransferFunction() {
        return new TransferFunction(this);
    }

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
    public Complex evaluateAt(double w) {
        Complex num = Polynomial.polyvalFromRoots(zeros, Complex.fromImaginary(w));
        Complex den = Polynomial.polyvalFromRoots(poles, Complex.fromImaginary(w));
        num.divideEquals(den);
        num.multiplyEquals(gain);
        return num;
    }
}
