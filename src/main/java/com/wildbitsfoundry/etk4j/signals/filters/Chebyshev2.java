package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;
import static com.wildbitsfoundry.etk4j.signals.filters.Filters.*;

public class Chebyshev2 extends AnalogFilter {
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
        double k = RationalFunction.calculateGain(zeros, poles);
        return new ZeroPoleGain(zeros, poles, k);
    }

    public static FilterOrderResults.OrderAndCutoffFrequency cheb2ord(FilterSpecs.LowPassSpecs specs) {
        return lowPassFilterOrder(specs, new Chebyshev2OrderCalculationStrategy());
    }

    public static FilterOrderResults.OrderAndCutoffFrequency cheb2ord(FilterSpecs.HighPassSpecs specs) {
        return highPassFilterOrder(specs, new Chebyshev2OrderCalculationStrategy());
    }

    public static FilterOrderResults.OrderAndCutoffFrequencies cheb2ord(FilterSpecs.BandPassSpecs specs) {
        return bandPassFilterOrder(specs, new Chebyshev2OrderCalculationStrategy());
    }

    public static FilterOrderResults.OrderAndCutoffFrequencies cheb2ord(FilterSpecs.BandStopSpecs specs) {
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
