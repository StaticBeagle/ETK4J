package com.wildbitsfoundry.etk4j.signals.filters;

import java.util.ArrayList;
import java.util.List;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;

import static com.wildbitsfoundry.etk4j.signals.filters.Filters.*;

public class Elliptic extends Filter {

    /**
     * Elliptic analog low pass filter prototype.
     * <br>
     * References: <pre>
     *     H. J. Orchard and A. N. Willson Jr., Elliptic functions for filter design,
     *     IEEE Trans. Circuits Syst., vol. 44, pp. 273–287, 1997
     * </pre>
     * @param n The order of the filter.
     * @param rp The pass band ripple in dB.
     * @param rs The stop band ripple in dB.
     * @return The zeros and poles of the elliptic filter.
     */
    public static ZeroPoleGain ellipAp(int n, double rp, double rs) {
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
            Complex c = mj.divide(new Complex(-u1, u2).cos());
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
        Complex num = ComplexArrays.product(zeros).multiply(Math.pow(-1, zeros.length));
        Complex den = ComplexArrays.product(poles).multiply(Math.pow(-1, poles.length));
        den.divideEquals(num);
        double k = den.real();
        if (n % 2 == 0) {
            double eps0 = e.get(0);
            k /= Math.sqrt(1 + eps0 * eps0);
        }
        return new ZeroPoleGain(zeros, poles, k);
    }

    /**
     * Elliptic filter order.
     * @param specs The filter design specifications.
     * @return The minimum order required to meet the design specifications.
     */
    public static LowPassResults ellipOrd(LowPassSpecs specs) {
        specs.validate();
        LowPassSpecs specsCopy = new LowPassSpecs(specs);
        double rp = 10 * Math.log10(specs.getPassBandRipple());
        double rs = 10 * Math.log10(specs.getStopBandAttenuation());
        specsCopy.setPassBandRipple(rp);
        specsCopy.setStopBandAttenuation(rs);
        return lowPassFilterOrder(specsCopy, new EllipticOrderCalculationStrategy());
    }

    /**
     * Elliptic filter order.
     * @param specs The filter design specifications.
     * @return The minimum order required to meet the design specifications.
     */
    public static HighPassResults ellipOrd(HighPassSpecs specs) {
        specs.validate();
        HighPassSpecs specsCopy = new HighPassSpecs(specs);
        double rp = 10 * Math.log10(specs.getPassBandRipple());
        double rs = 10 * Math.log10(specs.getStopBandAttenuation());
        specsCopy.setPassBandRipple(rp);
        specsCopy.setStopBandAttenuation(rs);
        return highPassFilterOrder(specsCopy, new EllipticOrderCalculationStrategy());
    }

    /**
     * Elliptic filter order.
     * @param specs The filter design specifications.
     * @return The minimum order required to meet the design specifications.
     */
    public static BandpassResults ellipOrd(BandpassSpecs specs) {
        specs.validate();
        BandpassSpecs specsCopy = new BandpassSpecs(specs);
        double rp = 10 * Math.log10(specs.getPassBandRipple());
        double rs = 10 * Math.log10(specs.getStopBandAttenuation());
        specsCopy.setPassBandRipple(rp);
        specsCopy.setStopBandAttenuation(rs);
        return bandpassFilterOrder(specsCopy, new EllipticOrderCalculationStrategy());
    }

    /**
     * Elliptic filter order.
     * @param specs The filter design specifications.
     * @return The minimum order required to meet the design specifications.
     */
    public static BandStopResults ellipOrd(BandStopSpecs specs) {
        specs.validate();
        BandStopSpecs specsCopy = new BandStopSpecs(specs);
        double rp = 10 * Math.log10(specs.getPassBandRipple());
        double rs = 10 * Math.log10(specs.getStopBandAttenuation());
        specsCopy.setPassBandRipple(rp);
        specsCopy.setStopBandAttenuation(rs);
        return bandStopFilterOrder(specsCopy, new EllipticOrderCalculationStrategy());
    }

    /**
     * Low pass filter realization.
     * @param n The order of the filter.
     * @param rp The pass band ripple.
     * @param rs The stop band ripple.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link Elliptic#newLowPassZPK(int, double, double, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newLowPass(int n, double rp, double rs, double wn) {
        ZeroPoleGain zpk = ellipAp(n, rp, rs);
        return lpToLp(zpk, wn);
    }

    /**
     * Low pass filter realization.
     * @param n The order of the filter.
     * @param rp The pass band ripple.
     * @param rs The stop band ripple.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newLowPassZPK(int n, double rp, double rs, double wn) {
        ZeroPoleGain zpk = ellipAp(n, rp, rs);
        return lpToLpZPK(zpk, wn);
    }

    /**
     * High pass filter realization.
     * @param n The order of the filter.
     * @param rp The pass band ripple.
     * @param rs The stop band ripple.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link Elliptic#newHighPassZPK(int, double, double, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newHighPass(int n, double rp, double rs, double wn) {
        ZeroPoleGain zpk = ellipAp(n, rp, rs);
        return lpToHp(zpk, wn);
    }

    /**
     * High pass filter realization.
     * @param n The order of the filter.
     * @param rp The pass band ripple.
     * @param rs The stop band ripple.
     * @param wn The cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newHighPassZPK(int n, double rp, double rs, double wn) {
        ZeroPoleGain zpk = ellipAp(n, rp, rs);
        return lpToHpZPK(zpk, wn);
    }

    /**
     * Bandpass filter realization.
     * @param n The order of the filter.
     * @param rp The pass band ripple.
     * @param rs The stop band ripple.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link Elliptic#newBandpassZPK(int, double, double, double, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newBandpass(int n, double rp, double rs, double wp1, double wp2) {
        ZeroPoleGain zpk = ellipAp(n, rp, rs);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBp(zpk, w0, bw);
    }

    /**
     * Bandpass filter realization.
     * @param n The order of the filter.
     * @param rp The pass band ripple.
     * @param rs The stop band ripple.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newBandpassZPK(int n, double rp, double rs, double wp1, double wp2) {
        ZeroPoleGain zpk = ellipAp(n, rp, rs);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBpZPK(zpk, w0, bw);
    }

    /**
     * Band stop filter realization.
     * @param n The order of the filter.
     * @param rp The pass band ripple.
     * @param rs The stop band ripple.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link TransferFunction} representation of the filter. {@link Elliptic#newBandStopZPK(int, double, double, double, double)},
     * is more numerically accurate than converting the zeros, poles, and gain into a {@link TransferFunction} so if the
     * coefficients of the numerator and denominator are not needed, please consider using the {@link ZeroPoleGain}
     * variant.
     */
    public static TransferFunction newBandStop(int n, double rp, double rs, double wp1, double wp2) {
        ZeroPoleGain zpk = ellipAp(n, rp, rs);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBs(zpk, w0, bw);
    }

    /**
     * Band stop filter realization.
     * @param n The order of the filter.
     * @param rp The pass band ripple.
     * @param rs The stop band ripple.
     * @param wp1 The lower cutoff frequency of the filter.
     * @param wp2 The upper cutoff frequency of the filter.
     * @return A {@link ZeroPoleGain} representation of the filter.
     */
    public static ZeroPoleGain newBandStopZPK(int n, double rp, double rs, double wp1, double wp2) {
        ZeroPoleGain zpk = ellipAp(n, rp, rs);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpToBsZPK(zpk, w0, bw);
    }

    private static void validateInputsLowPass(int n, double rp, double rs, double wn) {
        if(n <= 0) {
            throw new IllegalArgumentException("The filter order n must be greater than zero.");
        }
        if(rp <= 0) {
            throw new IllegalArgumentException("The pass band ripple rp must be greater than zero.");
        }
        if(rs <= 0) {
            throw new IllegalArgumentException("The stop band ripple rs must be greater than zero.");
        }
        if(wn <= 0) {
            throw new IllegalArgumentException("The cutoff frequency wn must be greater than zero.");
        }
    }

    private static void validateInputsHighPass(int n, double rp, double rs, double wn) {
        validateInputsLowPass(n, rp, rs, wn);
    }

    private static void validateInputsBandpass(int n, double rp, double rs, double wp1, double wp2) {
        if (n <= 0) {
            throw new IllegalArgumentException("The filter order n must be greater than zero.");
        }
        if(rp <= 0) {
            throw new IllegalArgumentException("The pass band ripple rp must be greater than zero.");
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

    private static void validateInputsBandStop(int n, double rp, double rs, double wp1, double wp2) {
        validateInputsBandpass(n, rp, rs, wp1, wp2);
    }
}
