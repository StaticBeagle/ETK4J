package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;
import com.wildbitsfoundry.etk4j.util.Tuples;

import static com.wildbitsfoundry.etk4j.signals.filters.AnalogFilter.lpTobp;

public class Chebyshev2 extends Filter {
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
                poles[i] = Complex.newComplex(-sinha * Math.cos(phik * pid), cosha * Math.sin(phik * pid));
                poles[i].divideEquals(Math.pow(poles[i].abs(), 2));
            }
        } else {
            for (int k = -(n - 1) >> 1, i = 0; k <= (n - 1) >> 1; ++k, ++i) {
                double phik = 180.0 * k * nInv;
                poles[i] = Complex.newComplex(-sinha * Math.cos(phik * pid), cosha * Math.sin(phik * pid));
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

    public static Tuples.Tuple2<Integer, Double> cheb2ord(FilterSpecs.LowPassSpecs specs) {
        return lowPassFilterOrder(specs, new Chebyshev2OrderCalculationStrategy());
    }

    public static Tuples.Tuple2<Integer, Double> cheb2ord(FilterSpecs.HighPassSpecs specs) {
        return highPassFilterOrder(specs, new Chebyshev2OrderCalculationStrategy());
    }

    public static Tuples.Tuple3<Integer, Double, Double> cheb2ord(FilterSpecs.BandPassSpecs specs) {
        return bandPassFilterOrder(specs, new Chebyshev2OrderCalculationStrategy());
    }

    public static Tuples.Tuple3<Integer, Double, Double> cheb2ord(FilterSpecs.BandStopSpecs specs) {
        return bandStopFilterOrder(specs, new Chebyshev2OrderCalculationStrategy());
    }


    public static TransferFunction newLowPass(int n, double rs, double wn) {
        ZeroPoleGain zpk = cheb2ap(n, rs);
        return AnalogFilter.lpTolp(zpk, wn);
    }

    public static TransferFunction newHighPass(int n, double rs, double wn) {
        ZeroPoleGain zpk = cheb2ap(n, rs);
        return AnalogFilter.lpTohp(zpk, wn);
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
        return AnalogFilter.lpTobs(zpk, w0, bw);
    }

    public static void main(String[] args) {
        double rs = 60;
        ZeroPoleGain zpk = cheb2ap(5, rs);
        FilterSpecs.LowPassSpecs lpSpecs = new FilterSpecs.LowPassSpecs();
        lpSpecs.setPassBandRipple(1.5); // 1.5 dB gain/ripple refer to note
        lpSpecs.setStopBandAttenuation(rs); // 60 dB at the stop band
        lpSpecs.setPassBandFrequency(2500); // 2500 Hz cutoff frequency
        lpSpecs.setStopBandFrequency(10000); // 10000 Hz stop band frequency
        Tuples.Tuple2<Integer, Double> cheb2ord = cheb2ord(lpSpecs);
        System.out.println(cheb2ord.getItem1());
        System.out.println(cheb2ord.getItem2());

        TransferFunction tf = newLowPass(cheb2ord.getItem1(), rs, cheb2ord.getItem2());
        System.out.println(tf);
        System.out.println();

        FilterSpecs.HighPassSpecs hpSpecs = new FilterSpecs.HighPassSpecs();
        hpSpecs.setPassBandRipple(0.2); // 0.2 dB gain/ripple refer to note
        hpSpecs.setStopBandAttenuation(rs); // 60 dB at the stop band
        hpSpecs.setPassBandFrequency(12); // 12 Hz cutoff frequency
        hpSpecs.setStopBandFrequency(0.2); // 0.2 Hz stop band frequency
        cheb2ord = cheb2ord(hpSpecs);
        System.out.println(cheb2ord.getItem1());
        System.out.println(cheb2ord.getItem2());

        tf = newHighPass(cheb2ord.getItem1(), rs, cheb2ord.getItem2());
        System.out.println(tf);
        System.out.println();

        FilterSpecs.BandPassSpecs bpSpecs = new FilterSpecs.BandPassSpecs();
        // The bandwidth of the filter starts at the LowerPassBandFrequency and
        // ends at the UpperPassBandFrequency. The filter has lower stop band
        // which is set LowerStopBandFrequency and the upper stop band can be set
        // with UpperStopBandFrequency. The attenuation at the stop bands can be
        // set with the LowerStopBandAttenuation and UpperStopBandAttenuation
        // respectively. In a frequency spectrum, the order of the frequencies will be:
        // LowerStopBandFrequency < LowerPassBandFrequency < UpperPassBandFrequency <
        // UpperStopBandFrequency
        bpSpecs.setLowerPassBandFrequency(190.0); // 190 Hz lower pass band frequency
        bpSpecs.setUpperPassBandFrequency(210.0); // 210 Hz upper pass band frequency
        bpSpecs.setLowerStopBandFrequency(180.0); // 180 Hz lower stop band frequency
        bpSpecs.setUpperStopBandFrequency(220.0); // 220 Hz upper stop band frequency
        bpSpecs.setPassBandRipple(0.2); // 0.2 dB gain/ripple refer to note
        bpSpecs.setLowerStopBandAttenuation(20.0); // 20 dB attenuation at the lower end of the skirt
        bpSpecs.setUpperStopBandAttenuation(20.0); // 20 dB attenuation at the upper end of the skirt
        bpSpecs.setStopBandAttenuation(20);
        Tuples.Tuple3<Integer, Double, Double> cheb2ordBP = cheb2ord(bpSpecs);

        tf = newBandPass(cheb2ordBP.getItem1(), 20, cheb2ordBP.getItem2(), cheb2ordBP.getItem3());
        System.out.println(tf);
        System.out.println();

        rs = 38.0;
        FilterSpecs.BandStopSpecs bsSpecs = new FilterSpecs.BandStopSpecs();
        // The notch of the filter starts at the LowerStopBandFrequency and
        // ends at the UpperStopBandFrequency. The filter has lower pass band
        // which is set LowerPassBandFrequency and the upper pass band can be set
        // with UpperPassBandFrequency. The attenuation at the notch can be
        // set with the StopBandAttenuation parameter and the attenuation/ripple
        // in the pass band can be set with the PassBandRipple parameter.
        // In a frequency spectrum, the order of the frequencies will be:
        // LowerPassBandFrequency < LowerStopBandFrequency < UpperStopBandFrequency <
        // UpperPassBandFrequency
        bsSpecs.setLowerPassBandFrequency(3.6e3); // 3600 Hz lower pass band frequency
        bsSpecs.setUpperPassBandFrequency(9.1e3); // 9100 Hz lower pass band frequency
        bsSpecs.setLowerStopBandFrequency(5.45e3); // 5450 Hz lower stop band frequency
        bsSpecs.setUpperStopBandFrequency(5.90e3); // 5900 Hz upper stop band frequency
        bsSpecs.setPassBandRipple(0.5); // 1.5 dB gain/ripple refer to note
        bsSpecs.setStopBandAttenuation(38.0); // 38 db attenuation at the notch
        Tuples.Tuple3<Integer, Double, Double> cheb2ordBS = cheb2ord(bsSpecs);

        tf = newBandStop(cheb2ordBS.getItem1(), 38.0, cheb2ordBS.getItem2(), cheb2ordBS.getItem3());
        System.out.println(tf);
        System.out.println();
    }
}
