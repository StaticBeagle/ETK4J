package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.LowPassSpecs;

public final class ButterWorth extends Filter {

    public static ZeroPoleGain buttAp(int n) {
        final double pid = Math.PI / 180.0;
        final double nInv = 1.0 / n;
        Complex[] poles = new Complex[n];
        if (n % 2 == 0) {
            for (int k = (-n >> 1) + 1, i = 0; k <= n >> 1; ++k, ++i) {
                double phik = nInv * (180.0 * k - 90.0);
                poles[n - i - 1] = Complex.newComplex(-Math.cos(phik * pid), Math.sin(phik * pid));
            }
        } else {
            for (int k = -(n - 1) >> 1, i = 0; k <= (n - 1) >> 1; ++k, ++i) {
                double phik = nInv * 180.0 * k;
                poles[n - i - 1] = Complex.newComplex(-Math.cos(phik * pid), Math.sin(phik * pid));
            }
        }
        Complex[] zeros = new Complex[0];
        double k = RationalFunction.calculateGain(zeros, poles);
        return new ZeroPoleGain(zeros, poles, k);
    }

    public static FilterOrderResults.OrderAndCutoffFrequency buttord(LowPassSpecs specs) {
        // TODO validate inputs
        return lowPassFilterOrder(specs, new ButterworthOrderCalculationStrategy());
    }

    public static FilterOrderResults.OrderAndCutoffFrequency buttord(HighPassSpecs specs) {
        return highPassFilterOrder(specs, new ButterworthOrderCalculationStrategy());
    }

    public static FilterOrderResults.OrderAndCutoffFrequencies buttord(BandPassSpecs specs) {
        return bandPassFilterOrder(specs, new ButterworthOrderCalculationStrategy());
    }

    public static FilterOrderResults.OrderAndCutoffFrequencies buttord(BandStopSpecs specs) {
        return bandStopFilterOrder(specs, new ButterworthOrderCalculationStrategy());
    }

    public static NumeratorDenominatorPair newLowPass(int n, double wn) {
        ZeroPoleGain zpk = buttAp(n);
        return AnalogFilter.lpTolp(zpk, wn);
    }

    public static NumeratorDenominatorPair newHighPass(int n, double wn) {
        ZeroPoleGain zpk = buttAp(n);
        return AnalogFilter.lpTohp(zpk, wn);
    }

    // TODO create exceptions. add checks to other filters
    public static NumeratorDenominatorPair newBandPass(int n, double wp1, double wp2) {
        if(n <= 0) {
            // throw
        }
        if(wp1 <= 0 || wp2 <= 0) {
           // throw
        }
        if(wp1 <= wp2) {
            // throw
        }
        ZeroPoleGain zpk = buttAp(n);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return AnalogFilter.lpTobp(zpk, w0, bw);
    }

    public static NumeratorDenominatorPair newBandStop(int n, double wp1, double wp2) {
        ZeroPoleGain zpk = buttAp(n);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return AnalogFilter.lpTobs(zpk, w0, bw);
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
