package com.wildbitsfoundry.etk4j.signals.filters;

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

    public static FilterOrderResults.OrderAndCutoffFrequency ellipord(FilterSpecs.LowPassSpecs specs) {
        FilterSpecs.LowPassSpecs specsCopy = new FilterSpecs.LowPassSpecs(specs);
        double rp = 10 * Math.log10(specs.getPassBandRipple());
        double rs = 10 * Math.log10(specs.getStopBandAttenuation());
        specsCopy.setPassBandRipple(rp);
        specsCopy.setStopBandAttenuation(rs);
        return lowPassFilterOrder(specsCopy, new EllipticOrderCalculationStrategy());
    }

    public static FilterOrderResults.OrderAndCutoffFrequency ellipord(FilterSpecs.HighPassSpecs specs) {
        FilterSpecs.HighPassSpecs specsCopy = new FilterSpecs.HighPassSpecs(specs);
        double rp = 10 * Math.log10(specs.getPassBandRipple());
        double rs = 10 * Math.log10(specs.getStopBandAttenuation());
        specsCopy.setPassBandRipple(rp);
        specsCopy.setStopBandAttenuation(rs);
        return highPassFilterOrder(specsCopy, new EllipticOrderCalculationStrategy());
    }

    public static FilterOrderResults.OrderAndCutoffFrequencies ellipord(FilterSpecs.BandPassSpecs specs) {
        FilterSpecs.BandPassSpecs specsCopy = new FilterSpecs.BandPassSpecs(specs);
        double rp = 10 * Math.log10(specs.getPassBandRipple());
        double rs = 10 * Math.log10(specs.getStopBandAttenuation());
        specsCopy.setPassBandRipple(rp);
        specsCopy.setStopBandAttenuation(rs);
        return bandPassFilterOrder(specsCopy, new EllipticOrderCalculationStrategy());
    }

    public static FilterOrderResults.OrderAndCutoffFrequencies ellipord(FilterSpecs.BandStopSpecs specs) {
        FilterSpecs.BandStopSpecs specsCopy = new FilterSpecs.BandStopSpecs(specs);
        double rp = 10 * Math.log10(specs.getPassBandRipple());
        double rs = 10 * Math.log10(specs.getStopBandAttenuation());
        specsCopy.setPassBandRipple(rp);
        specsCopy.setStopBandAttenuation(rs);
        return bandStopFilterOrder(specsCopy, new EllipticOrderCalculationStrategy());
    }

    public static NumeratorDenominatorPair newLowPass(int n, double rp, double rs, double wn) {
        ZeroPoleGain zpk = ellipap(n, rp, rs);
        return AnalogFilter.lpTolp(zpk, wn);
    }

    public static NumeratorDenominatorPair newHighPass(int n, double rp, double rs, double wn) {
        ZeroPoleGain zpk = ellipap(n, rp, rs);
        return AnalogFilter.lpTohp(zpk, wn);
    }

    public static NumeratorDenominatorPair newBandPass(int n, double rp, double rs, double wp1, double wp2) {
        ZeroPoleGain zpk = ellipap(n, rp, rs);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return AnalogFilter.lpTobp(zpk, w0, bw);
    }

    public static NumeratorDenominatorPair newBandStop(int n, double rp, double rs, double wp1, double wp2) {
        ZeroPoleGain zpk = ellipap(n, rp, rs);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return AnalogFilter.lpTobs(zpk, w0, bw);
    }
}
