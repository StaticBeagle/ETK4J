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

    public static LowPassResults buttord(LowPassSpecs specs) {
        specs.validate();
        return lowPassFilterOrder(specs, new ButterworthOrderCalculationStrategy());
    }

    public static HighPassResults buttord(HighPassSpecs specs) {
        specs.validate();
        return highPassFilterOrder(specs, new ButterworthOrderCalculationStrategy());
    }

    public static BandpassResults buttord(BandpassSpecs specs) {
        specs.validate();
        return bandpassFilterOrder(specs, new ButterworthOrderCalculationStrategy());
    }

    public static BandStopResults buttord(BandStopSpecs specs) {
        specs.validate();
        return bandStopFilterOrder(specs, new ButterworthOrderCalculationStrategy());
    }

    public static TransferFunction newLowPass(int n, double wn) {
        ZeroPoleGain zpk = buttAp(n);
        return lpTolp(zpk, wn);
    }

    public static ZeroPoleGain newLowPassZPK(int n, double wn) {
        ZeroPoleGain zpk = buttAp(n);
        return lpTolpZPK(zpk, wn);
    }

    public static TransferFunction newHighPass(int n, double wn) {
        ZeroPoleGain zpk = buttAp(n);
        return lpTohp(zpk, wn);
    }

    // TODO create exceptions. add checks to other filters
    public static TransferFunction newBandPass(int n, double wp1, double wp2) {
        if (n <= 0) {
            // throw
        }
        if (wp1 <= 0 || wp2 <= 0) {
            // throw
        }
        if (wp1 <= wp2) {
            // throw
        }
        ZeroPoleGain zpk = buttAp(n);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpTobp(zpk, w0, bw);
    }

    public static ZeroPoleGain newBandPassZPK(int n, double wp1, double wp2) {
        if (n <= 0) {
            // throw
        }
        if (wp1 <= 0 || wp2 <= 0) {
            // throw
        }
        if (wp1 <= wp2) {
            // throw
        }
        ZeroPoleGain zpk = buttAp(n);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpTobpZPK(zpk, w0, bw);
    }

    public static TransferFunction newBandStop(int n, double wp1, double wp2) {
        ZeroPoleGain zpk = buttAp(n);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpTobs(zpk, w0, bw);
    }

    public static ZeroPoleGain newBandStopZPK(int n, double wp1, double wp2) {
        ZeroPoleGain zpk = buttAp(n);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpTobsZPK(zpk, w0, bw);
    }

//    private static void checkInputsLowPassHighPass(int n, double wn) {
//        if (n <= 0) {
//            // throw
//        }
//        if (wn <= 0) {
//            // throw
//        }
//    }
}
