package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;
import com.wildbitsfoundry.etk4j.util.Tuples;

import java.util.ArrayList;
import java.util.List;

public class Elliptic extends Filter {

    // H. J. Orchard and A. N. Willson Jr., “Elliptic functions for filter design,”
    // IEEE Trans. Circuits Syst., vol. 44, pp. 273–287, 1997
    public static ZeroPoleGain ellipap(int n, double rp, double rs) {
        if (n == 1) {
            // filter becomes Chebyshev I
            Complex[] z = new Complex[0];
            Complex[] p = new Complex[1];
            p[0] = Complex.fromReal(-Math.sqrt(1.0 / (Math.pow(10.0, rp * 0.1) - 1.0)));
            double k = -p[0].real();
            return new ZeroPoleGain(z, p, k);
        }

        double dbn = Math.log(10.0) * 0.05;
        int n0 = (int) MathETK.rem(n, 2);
        int n3 = (n - n0) >> 1;
        double apn = dbn * rp;
        double asn = dbn * rs;

        List<Double> e = new ArrayList<>();
        e.add(Math.sqrt(2.0 * Math.exp(apn) * Math.sinh(apn)));

        List<Double> g = new ArrayList<>();
        g.add(e.get(0) / Math.sqrt(Math.exp(2 * asn) - 1));

        double v = g.get(0);
        int m2 = 0;
        while (v > 1.0e-150) {
            v = (v / (1.0 + Math.sqrt(1 - v * v)));
            v *= v;
            ++m2;
            g.add(v);
        }

        int m1 = 0;
        List<Double> ek = new ArrayList<>(m1);
        for (int i = 0; i < 10; ++i) {
            m1 = m2 + i;
            while (ek.size() <= m1) {
                ek.add(0.0);
            }
            ek.set(m1, 4.0 * Math.pow((g.get(m2) / 4.0), Math.pow(2.0, i) / n));
            if (ek.get(m1) < 1.0e-14) {
                break;
            }
        }

        for (int en = m1; en >= 1; --en) {
            ek.set(en - 1, 2.0 * Math.sqrt(ek.get(en)) / (1.0 + ek.get(en)));
        }

        double a = 0.0;
        for (int en = 1; en <= m2; ++en) {
            a = (1.0 + g.get(en)) * e.get(en - 1) * 0.5;
            e.add(a + Math.sqrt(a * a + g.get(en)));
        }

        double u2 = Math.log((1 + Math.sqrt(1 + Math.pow(e.get(m2), 2))) / e.get(m2)) / n;
        Complex[] zeros = new Complex[n % 2 != 0 ? n - 1 : n];
        Complex[] poles = new Complex[n];
        Complex j = Complex.fromImaginary(1.0);
        Complex mj = j.conj();
        for (int i = 0, m = zeros.length - 1; i < n3; ++i, m = m - 2) {
            double u1 = (2.0 * i + 1.0) * Math.PI / (2.0 * n);
            Complex c = mj.divide(Complex.newComplex(-u1, u2).cos());
            double d = 1.0 / Math.cos(u1);
            for (int en = m1; en >= 1; --en) {
                double k = ek.get(en);
                c = c.subtract(c.invert().multiply(k));
                c.divideEquals(1 + k);
                d = (d + k / d) / (1 + k);
            }
            Complex pole = c.invert();
            poles[m] = pole;
            poles[m - 1] = pole.conj();
            Complex zero = Complex.fromImaginary(d / ek.get(0));
            zeros[m] = zero;
            zeros[m - 1] = zero.conj();
        }
        if (n0 == 1) {
            a = 1.0 / Math.sinh(u2);
            for (int en = m1; en >= 1; --en) {
                double k = ek.get(en);
                a = (a - k / a) / (1 + k);
            }
            poles[n - 1] = Complex.fromReal(-1.0 / a);
        }
        double k = RationalFunction.calculateGain(zeros, poles);
        if (n % 2 == 0) {
            double eps0 = e.get(0);
            k /= Math.sqrt(1 + eps0 * eps0);
        }
        return new ZeroPoleGain(zeros, poles, k);
    }

    public static Tuples.Tuple2<Integer, Double> ellipord(FilterSpecs.LowPassSpecs specs) {
        double rp = 10 * Math.log10(specs.getPassBandRipple());
        double rs = 10 * Math.log10(specs.getStopBandAttenuation());
        specs.setPassBandRipple(rp);
        specs.setStopBandAttenuation(rs);
        return lowPassFilterOrder(specs, new EllipticOrderCalculationStrategy());
    }

    public static Tuples.Tuple2<Integer, Double> ellipord(FilterSpecs.HighPassSpecs specs) {
        double rp = 10 * Math.log10(specs.getPassBandRipple());
        double rs = 10 * Math.log10(specs.getStopBandAttenuation());
        specs.setPassBandRipple(rp);
        specs.setStopBandAttenuation(rs);
        return highPassFilterOrder(specs, new EllipticOrderCalculationStrategy());
    }

    public static Tuples.Tuple3<Integer, Double, Double> ellipord(FilterSpecs.BandPassSpecs specs) {
        double rp = 10 * Math.log10(specs.getPassBandRipple());
        double rs = 10 * Math.log10(specs.getStopBandAttenuation());
        specs.setPassBandRipple(rp);
        specs.setStopBandAttenuation(rs);
        return bandPassFilterOrder(specs, new EllipticOrderCalculationStrategy());
    }

    // TODO
    public static Tuples.Tuple3<Integer, Double, Double> ellipord(FilterSpecs.BandStopSpecs specs) {
        double rp = 10 * Math.log10(specs.getPassBandRipple());
        double rs = 10 * Math.log10(specs.getStopBandAttenuation());
        specs.setPassBandRipple(rp);
        specs.setStopBandAttenuation(rs);
        return bandStopFilterOrder(specs, new EllipticOrderCalculationStrategy());
    }

    public static TransferFunction newLowPass(int n, double rp, double rs, double wn) {
        ZeroPoleGain zpk = ellipap(n, rp, rs);
        return AnalogFilter.lpTolp(zpk, wn);
    }

    public static TransferFunction newHighPass(int n, double rp, double rs, double wn) {
        ZeroPoleGain zpk = ellipap(n, rp, rs);
        return AnalogFilter.lpTohp(zpk, wn);
    }

    // TODO
    public static TransferFunction newBandPass(int n, double rp, double rs, double wp1, double wp2) {
        ZeroPoleGain zpk = ellipap(n, rp, rs);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return AnalogFilter.lpTobp(zpk, w0, bw);
    }

    // TODO
    public static TransferFunction newBandStop(int n, double rp, double rs, double wp1, double wp2) {
        ZeroPoleGain zpk = ellipap(n, rp, rs);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return AnalogFilter.lpTobs(zpk, w0, bw);
    }

    public static void main(String[] args) {
        double rp = 1.5;
        ZeroPoleGain zpk = ellipap(5, rp, 60);
        FilterSpecs.LowPassSpecs lpSpecs = new FilterSpecs.LowPassSpecs();
        lpSpecs.setPassBandRipple(rp); // 1.5 dB gain/ripple refer to note
        lpSpecs.setStopBandAttenuation(60.0); // 60 dB at the stop band
        lpSpecs.setPassBandFrequency(2500); // 2500 Hz cutoff frequency
        lpSpecs.setStopBandFrequency(10000); // 10000 Hz stop band frequency
        Tuples.Tuple2<Integer, Double> ellipord = ellipord(lpSpecs);
        System.out.println(ellipord.getItem1());
        System.out.println(ellipord.getItem2());

        TransferFunction tf = newLowPass(ellipord.getItem1(), rp, 60, ellipord.getItem2());
        System.out.println(tf);
        System.out.println();

        rp = 0.2;
        FilterSpecs.HighPassSpecs hpSpecs = new FilterSpecs.HighPassSpecs();
        hpSpecs.setPassBandRipple(rp); // 0.2 dB gain/ripple refer to note
        hpSpecs.setStopBandAttenuation(60.0); // 60 dB at the stop band
        hpSpecs.setPassBandFrequency(12); // 12 Hz cutoff frequency
        hpSpecs.setStopBandFrequency(0.2); // 0.2 Hz stop band frequency
        ellipord = ellipord(hpSpecs);
        System.out.println(ellipord.getItem1());
        System.out.println(ellipord.getItem2());

        tf = newHighPass(ellipord.getItem1(), rp, 60, ellipord.getItem2());
        System.out.println(tf);
        System.out.println();

        rp = 0.2;
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
        bpSpecs.setPassBandRipple(rp); // 0.2 dB gain/ripple refer to note
        bpSpecs.setLowerStopBandAttenuation(20.0); // 20 dB attenuation at the lower end of the skirt
        bpSpecs.setUpperStopBandAttenuation(20.0); // 20 dB attenuation at the upper end of the skirt
        bpSpecs.setStopBandAttenuation(20);
        Tuples.Tuple3<Integer, Double, Double> ellipordBP = ellipord(bpSpecs);

        tf = newBandPass(ellipordBP.getItem1(), rp, 20, ellipordBP.getItem2(), ellipordBP.getItem3());
        System.out.println(tf);
        System.out.println();

        rp = 0.5;
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
        bsSpecs.setPassBandRipple(rp); // 1.5 dB gain/ripple refer to note
        bsSpecs.setStopBandAttenuation(38.0); // 38 db attenuation at the notch
        Tuples.Tuple3<Integer, Double, Double> ellipordBS = ellipord(bsSpecs);

        tf = newBandStop(ellipordBS.getItem1(), rp, 38, ellipordBS.getItem2(), ellipordBS.getItem3());
        System.out.println(tf);
        System.out.println();
    }
}
