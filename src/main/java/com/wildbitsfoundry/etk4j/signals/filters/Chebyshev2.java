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
    public static ZeroPoleGain cheb2ap(int n, double rs) {
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

    public static LowPassResults cheb2ord(LowPassSpecs specs) {
        return lowPassFilterOrder(specs, new Chebyshev2OrderCalculationStrategy());
    }

    public static HighPassResults cheb2ord(HighPassSpecs specs) {
        return highPassFilterOrder(specs, new Chebyshev2OrderCalculationStrategy());
    }

    public static BandpassResults cheb2ord(BandpassSpecs specs) {
        return bandpassFilterOrder(specs, new Chebyshev2OrderCalculationStrategy());
    }

    public static BandStopResults cheb2ord(BandStopSpecs specs) {
        return bandStopFilterOrder(specs, new Chebyshev2OrderCalculationStrategy());
    }

    public static TransferFunction newLowPass(int n, double rs, double wn) {
        ZeroPoleGain zpk = cheb2ap(n, rs);
        return lpTolp(zpk, wn);
    }

    public static TransferFunction newHighPass(int n, double rs, double wn) {
        ZeroPoleGain zpk = cheb2ap(n, rs);
        return lpTohp(zpk, wn);
    }

    public static TransferFunction newBandPass(int n, double rs, double wp1, double wp2) {
        ZeroPoleGain zpk = cheb2ap(n, rs);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpTobp(zpk, w0, bw);
    }

    public static TransferFunction newBandStop(int n, double rs, double wp1, double wp2) {
        ZeroPoleGain zpk = cheb2ap(n, rs);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpTobs(zpk, w0, bw);
    }
}
