package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;

import static com.wildbitsfoundry.etk4j.signals.filters.Filters.*;

public final class ButterWorth extends AnalogFilter {

    /**
     * Butterworth analog low pass filter prototype.
     * <br>
     * References:
     * <pre>
     *     Rolf Schaumann and Mac E. Van Valkenburg, "Design Of Analog Filters"
     * </pre>
     *
     * @param n The order of the filter.
     * @return The zeros and poles of the Butterworth filter.
     */
    public static ZeroPoleGain buttAp(int n) {
        final double pid = Math.PI / 180.0;
        final double nInv = 1.0 / n;
        Complex[] poles = new Complex[n];
        if (n % 2 == 0) {
            for (int k = (-n >> 1) + 1, i = 0; k <= n >> 1; ++k, ++i) {
                double phik = nInv * (180.0 * k - 90.0);
                poles[n - i - 1] = new Complex(-Math.cos(phik * pid), Math.sin(phik * pid));
            }
        } else {
            for (int k = -(n - 1) >> 1, i = 0; k <= (n - 1) >> 1; ++k, ++i) {
                double phik = nInv * 180.0 * k;
                poles[n - i - 1] = new Complex(-Math.cos(phik * pid), Math.sin(phik * pid));
            }
        }
        Complex[] zeros = new Complex[0];
        Complex num = ComplexArrays.product(zeros).multiply(Math.pow(-1, zeros.length));
        Complex den = ComplexArrays.product(poles).multiply(Math.pow(-1, poles.length));
        den.divideEquals(num);
        return new ZeroPoleGain(zeros, poles, den.real());
    }

    /**
     * Butterworth filter order.
     * @param specs The filter design specifications.
     * @return The minimum order required to meet the design specifications.
     */
    public static LowPassResults buttOrd(LowPassSpecs specs) {
        specs.validate();
        return lowPassFilterOrder(specs, new ButterworthOrderCalculationStrategy());
    }

    /**
     * Butterworth filter order.
     * @param specs The filter design specifications.
     * @return The minimum order required to meet the design specifications.
     */
    public static HighPassResults buttOrd(HighPassSpecs specs) {
        specs.validate();
        return highPassFilterOrder(specs, new ButterworthOrderCalculationStrategy());
    }

    /**
     * Butterworth filter order.
     * @param specs The filter design specifications.
     * @return The minimum order required to meet the design specifications.
     */
    public static BandpassResults buttOrd(BandpassSpecs specs) {
        specs.validate();
        return bandpassFilterOrder(specs, new ButterworthOrderCalculationStrategy());
    }

    /**
     * Butterworth filter order.
     * @param specs The filter design specifications.
     * @return The minimum order required to meet the design specifications.
     */
    public static BandStopResults buttOrd(BandStopSpecs specs) {
        specs.validate();
        return bandStopFilterOrder(specs, new ButterworthOrderCalculationStrategy());
    }

    /**
     * Low pass filter realization.
     * @param n The order of the filter.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link ButterWorth#newLowPassZPK(int, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newLowPass(int n, double wn) {
        validateInputsLowPass(n, wn);
        ZeroPoleGain zpk = buttAp(n);
        return lpToLp(zpk, wn);
    }

    /**
     * Low pass filter realization.
     * @param n The order of the filter.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newLowPassZPK(int n, double wn) {
        validateInputsLowPass(n, wn);
        ZeroPoleGain zpk = buttAp(n);
        return lpToLpZPK(zpk, wn);
    }

    /**
     * High pass filter realization.
     * @param n The order of the filter.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link ButterWorth#newHighPassZPK(int, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newHighPass(int n, double wn) {
        validateInputsHighPass(n, wn);
        ZeroPoleGain zpk = buttAp(n);
        return lpToHp(zpk, wn);
    }

    /**
     * High pass filter realization.
     * @param n The order of the filter.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newHighPassZPK(int n, double wn) {
        validateInputsHighPass(n, wn);
        ZeroPoleGain zpk = buttAp(n);
        return lpToHpZPK(zpk, wn);
    }

    /**
     * Bandpass filter realization.
     * @param n The order of the filter.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link ButterWorth#newBandpassZPK(int, double, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newBandpass(int n, double wp1, double wp2) {
        validateInputsBandpass(n, wp1, wp2);
        ZeroPoleGain zpk = buttAp(n);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBp(zpk, w0, bw);
    }

    /**
     * Bandpass filter realization.
     * @param n The order of the filter.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newBandpassZPK(int n, double wp1, double wp2) {
        validateInputsBandpass(n, wp1, wp2);
        ZeroPoleGain zpk = buttAp(n);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBpZPK(zpk, w0, bw);
    }

    /**
     * Band stop filter realization.
     * @param n The order of the filter.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link ButterWorth#newBandStopZPK(int, double, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newBandStop(int n, double wp1, double wp2) {
        validateInputsBandStop(n, wp1, wp2);
        ZeroPoleGain zpk = buttAp(n);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBs(zpk, w0, bw);
    }

    /**
     * Band stop filter realization.
     * @param n The order of the filter.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newBandStopZPK(int n, double wp1, double wp2) {
        validateInputsBandStop(n, wp1, wp2);
        ZeroPoleGain zpk = buttAp(n);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBsZPK(zpk, w0, bw);
    }

    protected static void validateInputsLowPass(int n, double wn) {
        if(n <= 0) {
            throw new IllegalArgumentException("The filter order n must be greater than zero.");
        }
        if(wn <= 0) {
            throw new IllegalArgumentException("The cutoff frequency wn must be greater than zero.");
        }
    }

    protected static void validateInputsHighPass(int n, double wn) {
        validateInputsLowPass(n, wn);
    }

    protected static void validateInputsBandpass(int n, double wp1, double wp2) {
        if (n <= 0) {
            throw new IllegalArgumentException("The filter order n must be greater than zero.");
        }
        if (wp1 <= 0 || wp2 <= 0) {
            throw new IllegalArgumentException("The cutoff frequencies wp1 & wp2 must be greater than zero.");
        }
        if (wp1 >= wp2) {
            throw new IllegalArgumentException("The cutoff frequency wp2 must be greater than wp1.");
        }
    }

    protected static void validateInputsBandStop(int n, double wp1, double wp2) {
        validateInputsBandpass(n, wp1, wp2);
    }
}
