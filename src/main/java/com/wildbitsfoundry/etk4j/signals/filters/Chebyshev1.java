package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;

import static com.wildbitsfoundry.etk4j.signals.filters.Filters.*;

public class Chebyshev1 extends AnalogFilter {

    public static ZeroPoleGain cheb1ap(int n, double rp) {
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

    public static FilterOrderResults.OrderAndCutoffFrequency cheb1ord(FilterSpecs.LowPassSpecs specs) {
        return lowPassFilterOrder(specs, new Chebyshev1OrderCalculationStrategy());
    }

    public static FilterOrderResults.OrderAndCutoffFrequency cheb1ord(FilterSpecs.HighPassSpecs specs) {
        return highPassFilterOrder(specs, new Chebyshev1OrderCalculationStrategy());
    }

    public static FilterOrderResults.OrderAndCutoffFrequencies cheb1ord(FilterSpecs.BandPassSpecs specs) {
        return bandPassFilterOrder(specs, new Chebyshev1OrderCalculationStrategy());
    }

    public static FilterOrderResults.OrderAndCutoffFrequencies cheb1ord(FilterSpecs.BandStopSpecs specs) {
        return bandStopFilterOrder(specs, new Chebyshev1OrderCalculationStrategy());
    }

    public static TransferFunction newLowPass(int n, double rp, double wn) {
        ZeroPoleGain zpk = cheb1ap(n, rp);
        return lpTolp(zpk, wn);
    }

    public static TransferFunction newHighPass(int n, double rp, double wn) {
        ZeroPoleGain zpk = cheb1ap(n, rp);
        return lpTohp(zpk, wn);
    }

    public static TransferFunction newBandPass(int n, double rp, double wp1, double wp2) {
        ZeroPoleGain zpk = cheb1ap(n, rp);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpTobp(zpk, w0, bw);
    }

    public static TransferFunction newBandStop(int n, double rp, double wp1, double wp2) {
        ZeroPoleGain zpk = cheb1ap(n, rp);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpTobs(zpk, w0, bw);
    }
}
