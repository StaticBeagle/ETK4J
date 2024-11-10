package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;

import static com.wildbitsfoundry.etk4j.signals.filters.Filters.*;

public class Chebyshev2 extends AnalogFilter {
    /**
     * Chebyshev type II (a.k.a. Inverse Chebyshev) analog low pass filter prototype.
     * <br>
     * References:
     * <pre>
     *     Rolf Schaumann and Mac E. Van Valkenburg, "Design Of Analog Filters"
     * </pre>
     *
     * @param n  The order of the filter.
     * @param rs The stop band ripple in dB.
     * @return The zeros and poles of the Chebyshev type II filter.
     */
    public static ZeroPoleGain cheb2Ap(int n, double rs) {
        double eps = 1.0 / Math.sqrt(Math.pow(10, rs * 0.1) - 1);

        double a = 1.0 / n * MathETK.asinh(1 / eps);
        double sinha = Math.sinh(a);
        double cosha = Math.cosh(a);

        Complex[] poles = new Complex[n];
        final double pid = Math.PI / 180.0;
        final double nInv = 1.0 / n;
        if (n % 2 == 0) {
            for (int k = (-n >> 1) + 1, i = 0; k <= n >> 1; ++k, ++i) {
                double phik = nInv * (180.0 * k - 90.0);
                poles[i] = new Complex(-sinha * Math.cos(phik * pid), cosha * Math.sin(phik * pid));
                poles[i].divideEquals(Math.pow(poles[i].abs(), 2));
            }
        } else {
            for (int k = -(n - 1) >> 1, i = 0; k <= (n - 1) >> 1; ++k, ++i) {
                double phik = 180.0 * k * nInv;
                poles[i] = new Complex(-sinha * Math.cos(phik * pid), cosha * Math.sin(phik * pid));
                poles[i].divideEquals(Math.pow(poles[i].abs(), 2));
            }
        }

        Complex[] zeros = new Complex[n % 2 == 0 ? n : n - 1];
        for (int k = 0; k < zeros.length; ) {
            Complex zero = Complex.fromImaginary(-1.0 / Math.cos(0.5 * Math.PI * (k + 1) * nInv));
            zeros[k++] = zero;
            zeros[k++] = zero.conj();
        }
        Complex num = ComplexArrays.product(zeros).multiply(Math.pow(-1, zeros.length));
        Complex den = ComplexArrays.product(poles).multiply(Math.pow(-1, poles.length));
        den.divideEquals(num);
        double k = den.real();
        return new ZeroPoleGain(zeros, poles, k);
    }

    /**
     * Inverse Chebyshev filter order.
     * @param specs The filter design specifications.
     * @return The minimum order required to meet the design specifications.
     */
    public static LowPassResults cheb2Ord(LowPassSpecs specs) {
        specs.validate();
        return lowPassFilterOrder(specs, new Chebyshev2OrderCalculationStrategy());
    }

    /**
     * Inverse Chebyshev filter order.
     * @param specs The filter design specifications.
     * @return The minimum order required to meet the design specifications.
     */
    public static HighPassResults cheb2Ord(HighPassSpecs specs) {
        specs.validate();
        return highPassFilterOrder(specs, new Chebyshev2OrderCalculationStrategy());
    }

    /**
     * Inverse Chebyshev filter order.
     * @param specs The filter design specifications.
     * @return The minimum order required to meet the design specifications.
     */
    public static BandpassResults cheb2Ord(BandpassSpecs specs) {
        specs.validate();
        return bandpassFilterOrder(specs, new Chebyshev2OrderCalculationStrategy());
    }

    /**
     * Inverse Chebyshev filter order.
     * @param specs The filter design specifications.
     * @return The minimum order required to meet the design specifications.
     */
    public static BandStopResults cheb2Ord(BandStopSpecs specs) {
        specs.validate();
        return bandStopFilterOrder(specs, new Chebyshev2OrderCalculationStrategy());
    }

    /**
     * Low pass filter realization.
     * @param n The order of the filter.
     * @param rs The stop band ripple.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link Chebyshev2#newLowPassZPK(int, double, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newLowPass(int n, double rs, double wn) {
        validateInputsLowPass(n, rs, wn);
        ZeroPoleGain zpk = cheb2Ap(n, rs);
        return lpToLp(zpk, wn);
    }

    /**
     * Low pass filter realization.
     * @param n The order of the filter.
     * @param rs The stop band ripple.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newLowPassZPK(int n, double rs, double wn) {
        validateInputsLowPass(n, rs, wn);
        ZeroPoleGain zpk = cheb2Ap(n, rs);
        return lpToLpZPK(zpk, wn);
    }

    /**
     * High pass filter realization.
     * @param n The order of the filter.
     * @param rs The stop band ripple.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link Chebyshev2#newHighPassZPK(int, double, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newHighPass(int n, double rs, double wn) {
        validateInputsHighPass(n, rs, wn);
        ZeroPoleGain zpk = cheb2Ap(n, rs);
        return lpToHp(zpk, wn);
    }

    /**
     * High pass filter realization.
     * @param n The order of the filter.
     * @param rs The stop band ripple.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newHighPassZPK(int n, double rs, double wn) {
        validateInputsHighPass(n, rs, wn);
        ZeroPoleGain zpk = cheb2Ap(n, rs);
        return lpToHpZPK(zpk, wn);
    }

    /**
     * Bandpass filter realization.
     * @param n The order of the filter.
     * @param rs The stop band ripple.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link Chebyshev2#newBandpassZPK(int, double, double, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newBandpass(int n, double rs, double wp1, double wp2) {
        validateInputsBandpass(n, rs, wp1, wp2);
        ZeroPoleGain zpk = cheb2Ap(n, rs);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBp(zpk, w0, bw);
    }

    /**
     * Bandpass filter realization.
     * @param n The order of the filter.
     * @param rs The stop band ripple.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newBandpassZPK(int n, double rs, double wp1, double wp2) {
        validateInputsBandpass(n, rs, wp1, wp2);
        ZeroPoleGain zpk = cheb2Ap(n, rs);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBpZPK(zpk, w0, bw);
    }

    /**
     * Band stop filter realization.
     * @param n The order of the filter.
     * @param rs The stop band ripple.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link Chebyshev2#newBandStopZPK(int, double, double, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newBandStop(int n, double rs, double wp1, double wp2) {
        validateInputsBandStop(n, rs, wp1, wp2);
        ZeroPoleGain zpk = cheb2Ap(n, rs);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBs(zpk, w0, bw);
    }

    /**
     * Band stop filter realization.
     * @param n The order of the filter.
     * @param rs The stop band ripple.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newBandStopZPK(int n, double rs, double wp1, double wp2) {
        validateInputsBandStop(n, rs, wp1, wp2);
        ZeroPoleGain zpk = cheb2Ap(n, rs);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBsZPK(zpk, w0, bw);
    }

    private static void validateInputsLowPass(int n, double rs, double wn) {
        if(n <= 0) {
            throw new IllegalArgumentException("The filter order n must be greater than zero.");
        }
        if(rs <= 0) {
            throw new IllegalArgumentException("The stop band ripple rs must be greater than zero.");
        }
        if(wn <= 0) {
            throw new IllegalArgumentException("The cutoff frequency wn must be greater than zero.");
        }
    }

    private static void validateInputsHighPass(int n, double rs, double wn) {
        validateInputsLowPass(n, rs, wn);
    }

    private static void validateInputsBandpass(int n, double rs, double wp1, double wp2) {
        if (n <= 0) {
            throw new IllegalArgumentException("The filter order n must be greater than zero.");
        }
        if(rs <= 0) {
            throw new IllegalArgumentException("The stop band ripple rs must be greater than zero.");
        }
        if (wp1 <= 0 || wp2 <= 0) {
            throw new IllegalArgumentException("The cutoff frequencies wp1 & wp2 must be greater than zero.");
        }
        if (wp1 >= wp2) {
            throw new IllegalArgumentException("The cutoff frequency wp2 must be greater than wp1.");
        }
    }

    private static void validateInputsBandStop(int n, double rs, double wp1, double wp2) {
        validateInputsBandpass(n, rs, wp1, wp2);
    }
}
