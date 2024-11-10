package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;

import static com.wildbitsfoundry.etk4j.signals.filters.Filters.*;

public class Chebyshev1 extends AnalogFilter {
    /**
     * Chebyshev type I analog low pass filter prototype.
     * <br>
     * References:
     * <pre>
     *     Rolf Schaumann and Mac E. Van Valkenburg, "Design Of Analog Filters"
     * </pre>
     *
     * @param n The order of the filter.
     * @param rp The pass band ripple in dB.
     * @return The zeros and poles of the Chebyshev type I filter.
     */
    public static ZeroPoleGain cheb1Ap(int n, double rp) {
        double eps = Math.sqrt(Math.pow(10, rp * 0.1) - 1);

        double a = 1.0 / n * MathETK.asinh(1 / eps);
        double sinha = Math.sinh(a);
        double cosha = Math.cosh(a);

        Complex[] poles = new Complex[n];
        final double pid = Math.PI / 180.0;
        final double nInv = 1.0 / n;
        if (n % 2 == 0) {
            for (int k = (-n >> 1) + 1, i = 0; k <= n >> 1; ++k, ++i) {
                double phik = nInv * (180.0 * k - 90.0);
                poles[n - i - 1] = new Complex(-sinha * Math.cos(phik * pid), cosha * Math.sin(phik * pid));
            }
        } else {
            for (int k = -(n - 1) >> 1, i = 0; k <= (n - 1) >> 1; ++k, ++i) {
                double phik = 180.0 * k * nInv;
                poles[n - i - 1] = new Complex(-sinha * Math.cos(phik * pid), cosha * Math.sin(phik * pid));
            }
        }
        Complex[] zeros = new Complex[0];
        Complex num = ComplexArrays.product(zeros).multiply(Math.pow(-1, zeros.length));
        Complex den = ComplexArrays.product(poles).multiply(Math.pow(-1, poles.length));
        den.divideEquals(num);
        double k = den.real();
        if (n % 2 == 0) {
            k /= Math.sqrt(1.0 + eps * eps);
        }
        return new ZeroPoleGain(zeros, poles, k);
    }

    /**
     * Chebyshev filter order.
     * @param specs The filter design specifications.
     * @return The minimum order required to meet the design specifications.
     */
    public static LowPassResults cheb1Ord(LowPassSpecs specs) {
        specs.validate();
        return lowPassFilterOrder(specs, new Chebyshev1OrderCalculationStrategy());
    }

    /**
     * Chebyshev filter order.
     * @param specs The filter design specifications.
     * @return The minimum order required to meet the design specifications.
     */
    public static HighPassResults cheb1Ord(HighPassSpecs specs) {
        specs.validate();
        return highPassFilterOrder(specs, new Chebyshev1OrderCalculationStrategy());
    }

    /**
     * Chebyshev filter order.
     * @param specs The filter design specifications.
     * @return The minimum order required to meet the design specifications.
     */
    public static BandpassResults cheb1Ord(BandpassSpecs specs) {
        specs.validate();
        return bandpassFilterOrder(specs, new Chebyshev1OrderCalculationStrategy());
    }

    /**
     * Chebyshev filter order.
     * @param specs The filter design specifications.
     * @return The minimum order required to meet the design specifications.
     */
    public static BandStopResults cheb1Ord(BandStopSpecs specs) {
        specs.validate();
        return bandStopFilterOrder(specs, new Chebyshev1OrderCalculationStrategy());
    }

    /**
     * Low pass filter realization.
     * @param n The order of the filter.
     * @param rp The pass band ripple.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link Chebyshev1#newLowPassZPK(int, double, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newLowPass(int n, double rp, double wn) {
        validateInputsLowPass(n, rp, wn);
        ZeroPoleGain zpk = cheb1Ap(n, rp);
        return lpToLp(zpk, wn);
    }

    /**
     * Low pass filter realization.
     * @param n The order of the filter.
     * @param rp The pass band ripple.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newLowPassZPK(int n, double rp, double wn) {
        validateInputsLowPass(n, rp, wn);
        ZeroPoleGain zpk = cheb1Ap(n, rp);
        return lpToLpZPK(zpk, wn);
    }

    /**
     * High pass filter realization.
     * @param n The order of the filter.
     * @param rp The pass band ripple.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link Chebyshev1#newHighPassZPK(int, double, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newHighPass(int n, double rp, double wn) {
        validateInputsHighPass(n, rp, wn);
        ZeroPoleGain zpk = cheb1Ap(n, rp);
        return lpToHp(zpk, wn);
    }

    /**
     * High pass filter realization.
     * @param n The order of the filter.
     * @param rp The pass band ripple.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newHighPassZPK(int n, double rp, double wn) {
        validateInputsHighPass(n, rp, wn);
        ZeroPoleGain zpk = cheb1Ap(n, rp);
        return lpToHpZPK(zpk, wn);
    }

    /**
     * Bandpass filter realization.
     * @param n The order of the filter.
     * @param rp The pass band ripple.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link Chebyshev1#newBandpassZPK(int, double, double, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newBandpass(int n, double rp, double wp1, double wp2) {
        validateInputsBandpass(n, rp, wp1, wp2);
        ZeroPoleGain zpk = cheb1Ap(n, rp);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBp(zpk, w0, bw);
    }

    /**
     * Bandpass filter realization.
     * @param n The order of the filter.
     * @param rp The pass band ripple.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newBandpassZPK(int n, double rp, double wp1, double wp2) {
        validateInputsBandpass(n, rp, wp1, wp2);
        ZeroPoleGain zpk = cheb1Ap(n, rp);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBpZPK(zpk, w0, bw);
    }

    /**
     * Band stop filter realization.
     * @param n The order of the filter.
     * @param rp The pass band ripple.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link Chebyshev1#newBandStopZPK(int, double, double, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newBandStop(int n, double rp, double wp1, double wp2) {
        validateInputsBandStop(n, rp, wp1, wp2);
        ZeroPoleGain zpk = cheb1Ap(n, rp);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBs(zpk, w0, bw);
    }

    /**
     * Band stop filter realization.
     * @param n The order of the filter.
     * @param rp The pass band ripple.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newBandStopZPK(int n, double rp, double wp1, double wp2) {
        validateInputsBandStop(n, rp, wp1, wp2);
        ZeroPoleGain zpk = cheb1Ap(n, rp);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBsZPK(zpk, w0, bw);
    }

    private static void validateInputsLowPass(int n, double rp, double wn) {
        if(n <= 0) {
            throw new IllegalArgumentException("The filter order n must be greater than zero.");
        }
        if(rp <= 0) {
            throw new IllegalArgumentException("The pass band ripple rp must be greater than zero.");
        }
        if(wn <= 0) {
            throw new IllegalArgumentException("The cutoff frequency wn must be greater than zero.");
        }
    }

    private static void validateInputsHighPass(int n, double rp, double wn) {
        validateInputsLowPass(n, rp, wn);
    }

    private static void validateInputsBandpass(int n, double rp, double wp1, double wp2) {
        if (n <= 0) {
            throw new IllegalArgumentException("The filter order n must be greater than zero.");
        }
        if(rp <= 0) {
            throw new IllegalArgumentException("The pass band ripple rp must be greater than zero.");
        }
        if (wp1 <= 0 || wp2 <= 0) {
            throw new IllegalArgumentException("The cutoff frequencies wp1 & wp2 must be greater than zero.");
        }
        if (wp1 >= wp2) {
            throw new IllegalArgumentException("The cutoff frequency wp2 must be greater than wp1.");
        }
    }

    private static void validateInputsBandStop(int n, double rp, double wp1, double wp2) {
        validateInputsBandpass(n, rp, wp1, wp2);
    }
}
