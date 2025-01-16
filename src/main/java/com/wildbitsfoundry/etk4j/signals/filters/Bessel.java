package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import com.wildbitsfoundry.etk4j.math.optimize.solvers.NewtonRaphson;
import com.wildbitsfoundry.etk4j.math.optimize.solvers.NewtonRaphsonComplex;
import com.wildbitsfoundry.etk4j.math.optimize.solvers.SolverResults;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;
import com.wildbitsfoundry.etk4j.util.DoubleArrays;

import static com.wildbitsfoundry.etk4j.signals.filters.Filters.*;

import java.util.Arrays;
import java.util.function.Function;

public class Bessel extends Filter {
    /**
     * Frequency normalization for Bessel filters.
     * <pre>
     *         PHASE
     *             The filter is normalized such that the phase response reaches its
     *             midpoint at an angular (e.g., rad/s) cutoff frequency of 1. This
     *             happens for both low-pass and high-pass filters, so this is the
     *             "phase-matched" case. [6]
     *             The magnitude response asymptotes are the same as a Butterworth
     *             filter of the same order with a cutoff of `Wn`.
     *             This is the default, and matches MATLAB's implementation.
     *         DELAY
     *             The filter is normalized such that the group delay in the passband
     *             is 1 (e.g., 1 second). This is the "natural" type obtained by
     *             solving Bessel polynomials
     *         MAGNITUDE
     *             The filter is normalized such that the gain magnitude is -3 dB at
     *             angular frequency 1. This is called "frequency normalization" by
     *             Bond. [1]
     *      </pre>
     */
    public enum FrequencyNormalization {
        PHASE,
        DELAY,
        MAGNITUDE
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
     */

    /**
     * Bessel analog low pass filter prototype.This method is an alias for {@link Bessel#besselApPhaseNormalized(int)}.
     * The filter is normalized such that the phase response reaches its
     * midpoint at an angular (e.g., rad/s) cutoff frequency of 1. This
     * happens for both low-pass and high-pass filters, so this is the
     * "phase-matched" case. [6]
     * The magnitude response asymptotes are the same as a Butterworth
     * filter of the same order with a cutoff of `Wn`.
     * This is the default, and matches MATLAB's implementation.
     * Port from scipy.
     *
     * @param n The order of the filter.
     * @return {@link ZeroPoleGain} representation of the Bessel filter Analog prototype.
     */
    public static ZeroPoleGain besselAp(int n) {
        return besselApPhaseNormalized(n);
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
     */

    /**
     * Bessel analog low pass filter prototype.
     * The filter is normalized such that the phase response reaches its
     * midpoint at an angular (e.g., rad/s) cutoff frequency of 1. This
     * happens for both low-pass and high-pass filters, so this is the
     * "phase-matched" case. [6]
     * The magnitude response asymptotes are the same as a Butterworth
     * filter of the same order with a cutoff of `Wn`.
     * This is the default, and matches MATLAB's implementation.
     * Port from scipy.
     *
     * @param n The order of the filter.
     * @return {@link ZeroPoleGain} representation of the Bessel filter Analog prototype.
     */
    public static ZeroPoleGain besselApPhaseNormalized(int n) {
        if (n == 0) {
            return new ZeroPoleGain(new Complex[]{}, new Complex[]{}, 1.0);
        }
        Complex[] zeros = {};
        Complex[] poles = Arrays.stream(besselZeros(n)).map(Complex::invert).toArray(Complex[]::new);
        double k = 1.0;

        double aLast = Math.floor(fallingFactorial(2 * n, n) / Math.pow(2, n));
        for (int i = 0; i < n; ++i) {
            poles[i].multiplyEquals(Math.pow(10, -Math.log10(aLast) / n));
        }
        return new ZeroPoleGain(zeros, poles, k);
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
     */

    /**
     * Bessel analog low pass filter prototype.
     * The filter is normalized such that the group delay in the passband
     * is 1 (e.g., 1 second). This is the "natural" type obtained by
     * solving Bessel polynomials.
     *
     * @param n The order of the filter.
     * @return {@link ZeroPoleGain} representation of the Bessel filter Analog prototype.
     */
    public static ZeroPoleGain besselApDelayNormalized(int n) {
        if (n == 0) {
            return new ZeroPoleGain(new Complex[]{}, new Complex[]{}, 1.0);
        }
        Complex[] zeros = {};
        Complex[] poles = Arrays.stream(besselZeros(n)).map(Complex::invert).toArray(Complex[]::new);
        double aLast = Math.floor(fallingFactorial(2 * n, n) / Math.pow(2, n));
        return new ZeroPoleGain(zeros, poles, aLast);
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
     */

    /**
     * Bessel analog low pass filter prototype.
     * The filter is normalized such that the gain magnitude is -3 dB at
     * angular frequency 1. This is called "frequency normalization" by
     * Bond. [1]
     *
     * @param n The order of the filter.
     * @return {@link ZeroPoleGain} representation of the Bessel filter Analog prototype.
     */
    public static ZeroPoleGain besselApMagnitudeNormalized(int n) {
        if (n == 0) {
            return new ZeroPoleGain(new Complex[]{}, new Complex[]{}, 1.0);
        }
        Complex[] zeros = {};
        Complex[] poles = Arrays.stream(besselZeros(n)).map(Complex::invert).toArray(Complex[]::new);

        double aLast = Math.floor(fallingFactorial(2 * n, n) / Math.pow(2, n));
        double k = aLast;
        // 1.0 / Math.sqrt(2.0) = -3 db which is equal to 10 ^ (-3.0 / 20.0)
        double normFactor = normFactor(poles, k, 1.0 / Math.sqrt(2.0));
        for (int i = 0; i < n; ++i) {
            poles[i].divideEquals(normFactor);
        }
        k = Math.pow(normFactor, -n) * aLast;
        return new ZeroPoleGain(zeros, poles, k);
    }

    /**
     * Low pass filter realization.
     *
     * @param n  The order of the filter.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link Bessel#newLowPassZPK(int, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newLowPass(int n, double wn) {
        ButterWorth.validateInputsLowPass(n, wn);
        ZeroPoleGain zpk = besselAp(n);
        return lpToLp(zpk, wn);
    }

    /**
     * Low pass filter realization.
     *
     * @param n  The order of the filter.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newLowPassZPK(int n, double wn) {
        ButterWorth.validateInputsLowPass(n, wn);
        return besselAp(n);
    }

    /**
     * High pass filter realization.
     *
     * @param n  The order of the filter.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link Bessel#newHighPassZPK(int, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newHighPass(int n, double wn) {
        ButterWorth.validateInputsHighPass(n, wn);
        ZeroPoleGain zpk = besselAp(n);
        return lpToHp(zpk, wn);
    }

    /**
     * High pass filter realization.
     *
     * @param n  The order of the filter.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newHighPassZPK(int n, double wn) {
        ButterWorth.validateInputsHighPass(n, wn);
        ZeroPoleGain zpk = besselAp(n);
        return lpToHpZPK(zpk, wn);
    }

    /**
     * Bandpass filter realization.
     *
     * @param n   The order of the filter.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link Bessel#newBandpassZPK(int, double, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newBandpass(int n, double wp1, double wp2) {
        ButterWorth.validateInputsBandpass(n, wp1, wp2);
        ZeroPoleGain zpk = besselAp(n);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBp(zpk, w0, bw);
    }

    /**
     * Bandpass filter realization.
     *
     * @param n   The order of the filter.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newBandpassZPK(int n, double wp1, double wp2) {
        ButterWorth.validateInputsBandpass(n, wp1, wp2);
        ZeroPoleGain zpk = besselAp(n);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBpZPK(zpk, w0, bw);
    }

    /**
     * Band stop filter realization.
     *
     * @param n   The order of the filter.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link Bessel#newBandStopZPK(int, double, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newBandStop(int n, double wp1, double wp2) {
        ButterWorth.validateInputsBandStop(n, wp1, wp2);
        ZeroPoleGain zpk = besselAp(n);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBs(zpk, w0, bw);
    }

    /**
     * Band stop filter realization.
     *
     * @param n   The order of the filter.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newBandStopZPK(int n, double wp1, double wp2) {
        ButterWorth.validateInputsBandStop(n, wp1, wp2);
        ZeroPoleGain zpk = besselAp(n);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBsZPK(zpk, w0, bw);
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy
     */
    private static Complex[] besselZeros(int n) {
        if (n == 0) {
            return new Complex[0];
        }
        Complex[] x0 = camposZeros(n);
        Complex[] x = aberth(fx, fp, x0, 1e-15, 50);

        for (int i = 0; i < n; ++i) {
            SolverResults<Complex> sr = new NewtonRaphsonComplex(
                    z -> {
                        Complex[] result = new Complex[1];
                        com.wildbitsfoundry.etk4j.math.specialfunctions.
                                Bessel.nonexpbesska01(n + 0.5, z.invert(), result, new Complex[1]);
                        return result[0];
                    }, x[i])
                    .derivative(
                            z -> {
                                Complex[] a0 = new Complex[1];
                                com.wildbitsfoundry.etk4j.math.specialfunctions.
                                        Bessel.nonexpbesska01(n - 0.5, z.invert(), a0, new Complex[1]);
                                a0[0] = a0[0].divide(z.pow(2).multiply(2));

                                Complex[] a1 = new Complex[1];
                                com.wildbitsfoundry.etk4j.math.specialfunctions.
                                        Bessel.nonexpbesska01(n + 0.5, z.invert(), a1, new Complex[1]);
                                a1[0] = a1[0].divide(z.pow(2));

                                Complex[] a2 = new Complex[1];
                                com.wildbitsfoundry.etk4j.math.specialfunctions.
                                        Bessel.nonexpbesska01(n + 1.5, z.invert(), a2, new Complex[1]);
                                a2[0] = a2[0].divide(z.pow(2).multiply(2));

                                return a0[0].subtract(a1[0]).add(a2[0]);
                            })
                    .absTolerance(1e-15)
                    .relTolerance(0.0)
                    .iterationLimit(50)
                    .solve();
            x[i] = sr.getValue();
        }

        Complex[][] mean = new Complex[2][n];
        mean[0] = ComplexArrays.deepCopy(x);
        for (int i = 0; i < n; ++i) {
            mean[1][n - i - 1] = x[i].conj();
        }

        for (int j = 0; j < n; ++j) {
            Complex[] temp = new Complex[mean.length];
            for (int i = 0; i < mean.length; ++i) {
                temp[i] = mean[i][j];
            }
            x[j] = ComplexArrays.mean(temp);
        }

        // zeros should sum to -1
        if (ComplexArrays.sum(x).add(1.0).abs() > 1e-15) {
            throw new RuntimeException("Generated zeros are inaccurate.");
        }
        return x;
    }

    private static Function<Complex[], Complex[]> fx = x -> {
        final int length = x.length;
        Complex[][] result = new Complex[length][1];
        for (int i = 0; i < length; ++i) {
            com.wildbitsfoundry.etk4j.math.specialfunctions.
                    Bessel.nonexpbesska01(length + 0.5, x[i].invert(), result[i], new Complex[1]);
        }
        return ComplexArrays.flatten(result);
    };

    private static Function<Complex[], Complex[]> fp = x -> {
        final int length = x.length;
        Complex[][] result = new Complex[length][1];
        for (int i = 0; i < length; ++i) {
            Complex[] a0 = new Complex[1];
            com.wildbitsfoundry.etk4j.math.specialfunctions.
                    Bessel.nonexpbesska01(length - 0.5, x[i].invert(), a0, new Complex[1]);
            a0[0] = a0[0].divide(x[i].pow(2).multiply(2));

            Complex[] a1 = new Complex[1];
            com.wildbitsfoundry.etk4j.math.specialfunctions.
                    Bessel.nonexpbesska01(length + 0.5, x[i].invert(), a1, new Complex[1]);
            a1[0] = a1[0].divide(x[i].pow(2));

            Complex[] a2 = new Complex[1];
            com.wildbitsfoundry.etk4j.math.specialfunctions.
                    Bessel.nonexpbesska01(length + 1.5, x[i].invert(), a2, new Complex[1]);
            a2[0] = a2[0].divide(x[i].pow(2).multiply(2));

            result[i] = new Complex[]{a0[0].subtract(a1[0]).add(a2[0])};
        }
        return ComplexArrays.flatten(result);
    };

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy
     */
    private static Complex[] camposZeros(int n) {
        if (n == 1) {
            return new Complex[]{Complex.fromReal(-1.0)};
        }

        double s = Polynomial.polyVal(new double[]{1, -3, 0, 2, 0, 0}, n);
        double b3 = Polynomial.polyVal(new double[]{-8, 16}, n) / s;
        double b2 = Polynomial.polyVal(new double[]{12, -12, -24}, n) / s;
        double b1 = Polynomial.polyVal(new double[]{-2, -12, 24, 8}, n) / s;
        double b0 = Polynomial.polyVal(new double[]{-1, 5, 0, -6, 0}, n) / s;

        double r = Polynomial.polyVal(new double[]{1, 2, 0, 0}, n);
        double a1 = Polynomial.polyVal(new double[]{-6, -6}, n) / r;
        double a2 = 6 / r;

        double[] k = DoubleArrays.linSteps(1, n);
        double[] x = Polynomial.polyVal(new double[]{a2, a1, 0}, k);
        double[] y = Polynomial.polyVal(new double[]{b3, b2, b1, b0}, k);
        return ComplexArrays.zip(x, y);
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy
     */
    private static Complex[] aberth(Function<Complex[], Complex[]> fx, Function<Complex[], Complex[]> fp,
                                    Complex[] x0, double tol, int maxIter) {
        final int n = x0.length;
        Complex[] x = ComplexArrays.deepCopy(x0);

        Complex[] beta = new Complex[n];
        int iter = 0;
        while (iter < maxIter) {
            Complex[] num = fx.apply(x);
            Complex[] den = fp.apply(x);
            Complex[] alpha = new Complex[n];
            for (int i = 0; i < n; ++i) {
                alpha[i] = num[i].divide(den[i]).uminus();
            }

            for (int i = 0; i < n; ++i) {
                beta[i] = new Complex();
                for (int j = i + 1; j < n; ++j) {
                    beta[i].addEquals(x[i].subtract(x[j]).invert());
                }
                for (int j = 0; j < i; ++j) {
                    beta[i].addEquals(x[i].subtract(x[j]).invert());
                }
            }

            for (int i = 0; i < n; ++i) {
                x[i].addEquals(alpha[i].divide(alpha[i].multiply(beta[i]).add(1.0)));
            }

            if (!Arrays.stream(x).allMatch(Complex::isFinite)) {
                throw new NonFiniteRootsException("Non finite roots in Aberth's method.");
            }

            boolean done = true;
            for (int i = 0; i < n; ++i) {
                done &= alpha[i].abs() <= tol;
            }
            if (done) {
                break;
            }
            if (++iter > maxIter) {
                throw new MaximumNumberOfIterationsExceededException("Aberth maximum number of iterations exceeded.");
            }
        }
        return x;
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy
     */
    private static double fallingFactorial(int x, int n) {
        double val = 1.0;
        for (int i = x - n + 1; i < x + 1; ++i) {
            val *= i;
        }
        return val;
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy
     */
    private static double normFactor(Complex[] poles, double k, double magDrop) {
        UnivariateFunction gw = w -> {
            Complex prod = Complex.fromReal(1.0);
            for (int i = 0; i < poles.length; ++i) {
                prod.multiplyEquals(Complex.fromImaginary(w).subtract(poles[i]));
            }
            return prod.invert().multiply(k).abs();
        };
        UnivariateFunction cutoff = w -> gw.evaluateAt(w) - magDrop;

        double result = new NewtonRaphson(cutoff, 1.5)
                .absTolerance(1.48e-8)
                .relTolerance(0.0)
                .iterationLimit(50)
                .solve()
                .getValue();
        return result;
    }
}
