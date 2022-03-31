package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;

import java.util.Arrays;

/**
 * The {@Filters} class is a utility class providing static methods for {@link AnalogFilter}.
 */
public final class Filters {
    private Filters() {
    }

    /**
     * Change the cutoff frequency of a low pass filter.
     *
     * @param num The numerator coefficients of the filter in descending order.
     * @param den The denominator coefficients of the filter in descending order.
     * @param w0  The new cutoff frequency.
     * @return The {@link TransferFunction} of the new low pass filter with new cutoff.
     */
    public static TransferFunction lpTolp(double[] num, double[] den, double w0) {
        TransferFunction tf = new TransferFunction(num, den);
        tf.substituteInPlace(1.0 / w0);
        tf.normalize();
        return tf;
    }

    /**
     * Change the cutoff frequency of a low pass filter.
     *
     * @param zpk The {@link ZeroPoleGain} representation of the filter.
     * @param w0  The new cutoff frequency.
     * @return The new {@link ZeroPoleGain} representation of the filter with new cutoff frequency.
     */
    public static ZeroPoleGain lpTolpZPK(ZeroPoleGain zpk, double w0) {
        Complex[] zeros = zpk.getZeros();
        Complex[] poles = zpk.getPoles();
        double k = zpk.getGain();
        int degree = getRelativeDegree(zeros, poles);
        ComplexArrays.multiplyElementWiseInPlace(zeros, w0);
        ComplexArrays.multiplyElementWiseInPlace(poles, w0);
        k *= Math.pow(w0, degree);
        return new ZeroPoleGain(zeros, poles, k);
    }

    /**
     * Change the cutoff frequency of a low pass filter.
     *
     * @param zpk The {@link ZeroPoleGain} representation of the filter.
     * @param w0  The new cutoff frequency.
     * @return The {@link TransferFunction} of the new low pass filter with new cutoff.
     */
    public static TransferFunction lpTolp(ZeroPoleGain zpk, double w0) {
        return new TransferFunction(lpTolpZPK(zpk, w0));
    }

    /**
     * Transform a low pass filter into a high pass filter.
     *
     * @param num The numerator coefficients of the filter in descending order.
     * @param den The denominator coefficients of the filter in descending order.
     * @param w0  The cutoff frequency of the high pass filter.
     * @return The {@link TransferFunction} of the new high pass filter.
     */
    public static TransferFunction lpTohp(double[] num, double[] den, double w0) {
        final int numDegree = num.length - 1;
        final int denDegree = den.length - 1;
        final int filterOrder = Math.max(numDegree, denDegree) + 1;

        // Reverse coefficients then scale them by the
        // order of the denominator i.e. pad with zeros
        double[] hpNumerator = new double[filterOrder];
        for (int i = numDegree, j = 0; i >= 0; --i, ++j) {
            hpNumerator[j] = num[i] * Math.pow(w0, j);
        }

        // Reverse coefficients then scale them by the
        // order of the numerator i.e. pad with zeros
        double[] hpDenominator = new double[filterOrder];
        for (int i = denDegree, j = 0; i >= 0; --i, ++j) {
            hpDenominator[j] = den[i] * Math.pow(w0, j);
        }
        TransferFunction tf = new TransferFunction(hpNumerator, hpDenominator);
        tf.normalize();
        return tf;
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. see https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
    */

    /**
     * Transform a low pass filter into a high pass filter.
     *
     * @param zpk The {@link ZeroPoleGain} representation of the filter.
     * @param w0  The cutoff frequency of the high pass filter.
     * @return The new {@link ZeroPoleGain} representation of the new high pass filter.
     */
    public static ZeroPoleGain lpTohpZPK(ZeroPoleGain zpk, double w0) {
        Complex[] zeros = zpk.getZeros();
        Complex[] poles = zpk.getPoles();
        double k = zpk.getGain();
        int degree = getRelativeDegree(zeros, poles);

        Complex[] zhp = ComplexArrays.divideElementWise(w0, zeros);
        Complex[] php = ComplexArrays.divideElementWise(w0, poles);

        zhp = ComplexArrays.concatenate(zhp, ComplexArrays.zeros(degree));

        zeros = Arrays.stream(zeros).map(Complex::uminus).toArray(Complex[]::new);
        poles = Arrays.stream(poles).map(Complex::uminus).toArray(Complex[]::new);
        k *= ComplexArrays.product(zeros).divide(ComplexArrays.product(poles)).real();

        return new ZeroPoleGain(zhp, php, k);
    }

    /**
     * Transform a low pass filter into a high pass filter.
     *
     * @param zpk The {@link ZeroPoleGain} representation of the filter.
     * @param w0  The cutoff frequency of the high pass filter.
     * @return The {@link TransferFunction} of the new high pass filter.
     */
    public static TransferFunction lpTohp(ZeroPoleGain zpk, double w0) {
        return new TransferFunction(lpTohpZPK(zpk, w0));
    }

    /**
     * Transform a low pass filter into a bandpass filter.
     *
     * @param num The numerator coefficients of the filter in descending order.
     * @param den The denominator coefficients of the filter in descending order.
     * @param w0  The center frequency.
     * @param bw  The bandwidth of the filter.
     * @return The {@link TransferFunction} of the new bandpass filter.
     */
    public static TransferFunction lpTobp(double[] num, double[] den, double w0, double bw) {
        Polynomial s = new Polynomial(bw, 0.0);
        Polynomial s2w02 = new Polynomial(1.0, 0, w0 * w0);

        RationalFunction bp = new RationalFunction(num, den);
        bp.substituteInPlace(new RationalFunction(s2w02, s));

        TransferFunction tf = new TransferFunction(bp);
        tf.normalize();
        return tf;
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. see https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
    */

    /**
     * Transform a low pass filter into a bandpass filter.
     *
     * @param zpk The {@link ZeroPoleGain} representation of the filter.
     * @param w0  The center frequency of the filter.
     * @param bw  The bandwidth of the filter.
     * @return The {@link ZeroPoleGain} representation of the new bandpass filter.
     */
    public static ZeroPoleGain lpTobpZPK(ZeroPoleGain zpk, double w0, double bw) {
        Complex[] zeros = zpk.getZeros();
        Complex[] poles = zpk.getPoles();
        double k = zpk.getGain();
        Complex[] zlp = ComplexArrays.multiplyElementWise(zeros, bw * 0.5);
        Complex[] plp = ComplexArrays.multiplyElementWise(poles, bw * 0.5);

        int degree = getRelativeDegree(zeros, poles);

        Complex[] left = new Complex[zlp.length];
        Complex[] right = new Complex[zlp.length];
        for (int i = 0; i < zlp.length; ++i) {
            left[i] = zlp[i].pow(2.0).subtract(w0 * w0).sqrt().add(zlp[i]);
            right[i] = zlp[i].pow(2.0).subtract(w0 * w0).sqrt().uminus().add(zlp[i]);
            if (zlp[i].real() == 0.0) {
                left[i] = Complex.fromImaginary(left[i].imag());
                right[i] = Complex.fromImaginary(right[i].imag());
            }
        }
        Complex[] zbp = ComplexArrays.concatenate(left, right);

        left = new Complex[plp.length];
        right = new Complex[plp.length];
        for (int i = 0; i < plp.length; ++i) {
            left[i] = plp[i].pow(2.0).subtract(w0 * w0).sqrt().add(plp[i]);
            right[i] = plp[i].pow(2.0).subtract(w0 * w0).sqrt().uminus().add(plp[i]);
        }
        Complex[] pbp = ComplexArrays.concatenate(left, right);

        zbp = ComplexArrays.concatenate(zbp, ComplexArrays.zeros(degree));
        k *= Math.pow(bw, degree);

        return new ZeroPoleGain(zbp, pbp, k);
    }

    /**
     * Transform a low pass filter into a bandpass filter.
     *
     * @param zpk The {@link ZeroPoleGain} representation of the filter.
     * @param w0  The center frequency.
     * @param bw  The bandwidth.
     * @return The {@link TransferFunction} of the new bandpass filter.
     */
    public static TransferFunction lpTobp(ZeroPoleGain zpk, double w0, double bw) {
        return new TransferFunction(lpTobpZPK(zpk, w0, bw));
    }

    /**
     * Transform low pass filter into a band stop filter.
     *
     * @param num The numerator coefficients of the filter in descending order.
     * @param den The denominator coefficients of the filter in descending order.
     * @param w0  The center frequency.
     * @param bw  The bandwidth.
     * @return The {@link TransferFunction} of the new band stop filter.
     */
    public static TransferFunction lpTobs(double[] num, double[] den, double w0, double bw) {
        Polynomial s = new Polynomial(bw, 0.0);
        Polynomial s2w02 = new Polynomial(1.0, 0.0, w0 * w0);

        RationalFunction bp = new RationalFunction(num, den);
        bp.substituteInPlace(new RationalFunction(s, s2w02));

        return new TransferFunction(bp);
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. see https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
    */

    /**
     * Transform a low pass filter into a band stop filter.
     *
     * @param zpk The {@link ZeroPoleGain} representation of the filter.
     * @param w0  The center frequency.
     * @param bw  The bandwidth.
     * @return The {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain lpTobsZPK(ZeroPoleGain zpk, double w0, double bw) {
        Complex[] zeros = zpk.getZeros();
        Complex[] poles = zpk.getPoles();
        double k = zpk.getGain();

        int degree = getRelativeDegree(zeros, poles);

        Complex[] zhp = ComplexArrays.divideElementWise(bw * 0.5, zeros);
        Complex[] php = ComplexArrays.divideElementWise(bw * 0.5, poles);

        Complex[] left = new Complex[zhp.length];
        Complex[] right = new Complex[zhp.length];
        for (int i = 0; i < zhp.length; ++i) {
            left[i] = zhp[i].add(zhp[i].pow(2.0).subtract(w0 * w0).sqrt());
            right[i] = zhp[i].subtract(zhp[i].pow(2.0).subtract(w0 * w0).sqrt());
            if (zhp[i].real() == 0.0) {
                left[i] = Complex.fromImaginary(left[i].imag());
                right[i] = Complex.fromImaginary(right[i].imag());
            }
        }
        Complex[] zbs = ComplexArrays.concatenate(left, right);

        left = new Complex[php.length];
        right = new Complex[php.length];
        for (int i = 0; i < php.length; ++i) {
            left[i] = php[i].add(php[i].pow(2.0).subtract(w0 * w0).sqrt());
            right[i] = php[i].subtract(php[i].pow(2.0).subtract(w0 * w0).sqrt());
        }
        Complex[] pbs = ComplexArrays.concatenate(left, right);

        Complex[] full = new Complex[degree];
        Arrays.fill(full, Complex.fromImaginary(w0));
        zbs = ComplexArrays.concatenate(zbs, full);

        Arrays.fill(full, Complex.fromImaginary(-w0));
        zbs = ComplexArrays.concatenate(zbs, full);

        zeros = Arrays.stream(zeros).map(Complex::uminus).toArray(Complex[]::new);
        poles = Arrays.stream(poles).map(Complex::uminus).toArray(Complex[]::new);
        k *= ComplexArrays.product(zeros).divide(ComplexArrays.product(poles)).real();

        return new ZeroPoleGain(zbs, pbs, k);
    }

    /**
     * Transform a low pass filter into a band stop filter.
     * @param zpk The {@link ZeroPoleGain} representation of the filter.
     * @param w0 The center frequency.
     * @param bw The bandwidth.
     * @return The {@link TransferFunction} of the new band stop filter.
     */
    public static TransferFunction lpTobs(ZeroPoleGain zpk, double w0, double bw) {
        return new TransferFunction(lpTobsZPK(zpk, w0, bw));
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. see https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
    */
    private static int getRelativeDegree(Complex[] zeros, Complex[] poles) {
        int degree = poles.length - zeros.length;
        if (degree < 0) {
            throw new NegativeFilterOrderException("The number of poles for the filter is less than the number of zeros."
                    + " Please check your inputs.");
        }
        return degree;
    }
}
