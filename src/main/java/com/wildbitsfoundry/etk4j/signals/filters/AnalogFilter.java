package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.LowPassSpecs;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;

import java.util.Arrays;
import java.util.function.BiFunction;

import static com.wildbitsfoundry.etk4j.math.optimize.minimizers.Brent.brentsMinimizer;
import static com.wildbitsfoundry.etk4j.math.optimize.minimizers.GoldenSection.goldenSectionMinimizer;

public class AnalogFilter {

    private int _order;
    private TransferFunction _tf;

    private static int getRelativeDegree(Complex[] zeros, Complex[] poles) {
        int degree = poles.length - zeros.length;
        if(degree < 0) {
            // throw
        }
        return degree;
    }

    public static NumeratorDenominatorPair lpTolp(ZeroPoleGain zpk, double w0) {
        Complex[] zeros = zpk.getZeros();
        Complex[] poles = zpk.getPoles();
        double k = zpk.getGain();
        int degree = getRelativeDegree(zeros, poles);
        ComplexArrays.multiplyInPlace(zeros, w0);
        ComplexArrays.multiplyInPlace(poles, w0);
        k *= Math.pow(w0, degree);
        RationalFunction rf = new RationalFunction(zeros, poles, k);
        return new NumeratorDenominatorPair(rf.getNumerator().getCoefficients(), rf.getDenominator().getCoefficients());
    }

    public static NumeratorDenominatorPair lpTobp(ZeroPoleGain zpk, double w0, double bw) {
        Complex[] zeros = zpk.getZeros();
        Complex[] poles = zpk.getPoles();
        double k = zpk.getGain();
        Complex[] zlp = ComplexArrays.multiply(zeros, bw * 0.5);
        Complex[] plp = ComplexArrays.multiply(poles, bw * 0.5);

        // TODO make this a method
        int degree = poles.length - zeros.length;
        // TODO if degree < 0 throw

        Complex[] left = new Complex[zlp.length];
        Complex[] right = new Complex[zlp.length];
        for(int i = 0; i < zlp.length; ++i) {
            left[i] = zlp[i].pow(2.0).subtract(w0 * w0).sqrt().add(zlp[i]);
            right[i] = zlp[i].pow(2.0).subtract(w0 * w0).sqrt().uminus().add(zlp[i]);
            if(zlp[i].real() == 0.0) {
                left[i] = Complex.fromImaginary(left[i].imag());
                right[i] = Complex.fromImaginary(right[i].imag());
            }
        }
        Complex[] zbp = ComplexArrays.concat(left, right);

        left = new Complex[plp.length];
        right = new Complex[plp.length];
        for(int i = 0; i < plp.length; ++i) {
            left[i] = plp[i].pow(2.0).subtract(w0 * w0).sqrt().add(plp[i]);
            right[i] = plp[i].pow(2.0).subtract(w0 * w0).sqrt().uminus().add(plp[i]);
        }
        Complex[] pbp = ComplexArrays.concat(left, right);

        zbp = ComplexArrays.concat(zbp, ComplexArrays.zeros(degree));
        k *= Math.pow(bw, degree);

        RationalFunction rf = new RationalFunction(zbp, pbp, k);
        return new NumeratorDenominatorPair(rf.getNumerator().getCoefficients(), rf.getDenominator().getCoefficients());
    }

    public static NumeratorDenominatorPair lpTohp(ZeroPoleGain zpk, double w0) {
        Complex[] zeros = zpk.getZeros();
        Complex[] poles = zpk.getPoles();
        double k = zpk.getGain();
        int degree = poles.length - zeros.length;
        // TODO if degree < 0 throw

        Complex[] zhp = ComplexArrays.divide(w0, zeros);
        Complex[] php = ComplexArrays.divide(w0, poles);

        zhp = ComplexArrays.concat(zhp, ComplexArrays.zeros(degree));

        zeros = Arrays.stream(zeros).map(Complex::uminus).toArray(Complex[]::new);
        poles = Arrays.stream(poles).map(Complex::uminus).toArray(Complex[]::new);
        k *= ComplexArrays.product(zeros).divide(ComplexArrays.product(poles)).real();

        RationalFunction rf = new RationalFunction(zhp, php, k);
        return new NumeratorDenominatorPair(rf.getNumerator().getCoefficients(), rf.getDenominator().getCoefficients());
    }

    public static NumeratorDenominatorPair lpTobs(ZeroPoleGain zpk, double w0, double bw) {
        Complex[] zeros = zpk.getZeros();
        Complex[] poles = zpk.getPoles();
        double k = zpk.getGain();

        int degree = getRelativeDegree(zeros, poles);

        Complex[] zhp = ComplexArrays.divide(bw * 0.5, zeros);
        Complex[] php = ComplexArrays.divide(bw * 0.5, poles);

        Complex[] left = new Complex[zhp.length];
        Complex[] right = new Complex[zhp.length];
        for(int i = 0; i < zhp.length; ++i) {
            //left[i] = zhp[i].pow(2.0).subtract(w0 * w0).sqrt().add(zhp[i]);
            left[i] = zhp[i].add(zhp[i].pow(2.0).subtract(w0 * w0).sqrt());
            right[i] = zhp[i].subtract(zhp[i].pow(2.0).subtract(w0 * w0).sqrt());
            if(zhp[i].real() == 0.0) {
                left[i] = Complex.fromImaginary(left[i].imag());
                right[i] = Complex.fromImaginary(right[i].imag());
            }
        }
        Complex[] zbs = ComplexArrays.concat(left, right);

        left = new Complex[php.length];
        right = new Complex[php.length];
        for(int i = 0; i < php.length; ++i) {
            left[i] = php[i].add(php[i].pow(2.0).subtract(w0 * w0).sqrt());
            right[i] = php[i].subtract(php[i].pow(2.0).subtract(w0 * w0).sqrt());
        }
        Complex[] pbs = ComplexArrays.concat(left, right);

        Complex[] full = new Complex[degree];
        Arrays.fill(full, Complex.fromImaginary(w0));
        zbs = ComplexArrays.concat(zbs, full);

        Arrays.fill(full, Complex.fromImaginary(-w0));
        zbs = ComplexArrays.concat(zbs, full);

        zeros = Arrays.stream(zeros).map(Complex::uminus).toArray(Complex[]::new);
        poles = Arrays.stream(poles).map(Complex::uminus).toArray(Complex[]::new);
        k *= ComplexArrays.product(zeros).divide(ComplexArrays.product(poles)).real();

        RationalFunction rf = new RationalFunction(zbs, pbs, k);
        return new NumeratorDenominatorPair(rf.getNumerator().getCoefficients(), rf.getDenominator().getCoefficients());
    }

    static class LowPassPrototype {
        private TransferFunction _tf;

        public LowPassPrototype(ZeroPoleGain zpk) {
            _tf = zpkToTF(zpk);
        }
    }

    protected AnalogFilter(int order, TransferFunction tf) {
        _order = order;
        _tf = tf;
    }

    public double[] getNumerator() {
        return _tf.getNumerator().getCoefficients();
    }

    public double[] getDenominator() {
        return _tf.getDenominator().getCoefficients();
    }

    @Override
    public String toString() {
        return _tf.toString();
    }

    public int getOrder() {
        return _order;
    }

    public static AnalogFilter newLowPass(LowPassSpecs specs, ApproximationType type) {
        double fp = specs.getPassBandFrequency();
        double fs = specs.getStopBandFrequency();
        double ap = specs.getPassBandRipple();
        double as = specs.getStopBandAttenuation();

        double wp = 2 * Math.PI * fp;
        double ws = 2 * Math.PI * fs;
        final int n = type.getMinOrderNeeded(wp, ws, ap, as);
        LowPassPrototype lp = new LowPassPrototype(type.buildLowPassPrototype(n, ap, as));
        lp._tf = lpTolp(lp._tf.getNumerator(), lp._tf.getDenominator(), fp);
        lp._tf.normalize();
        return new AnalogFilter(n, lp._tf);
    }

    public static AnalogFilter newHighPass(HighPassSpecs specs, ApproximationType type) {
        double fp = specs.getPassBandFrequency();
        double fs = specs.getStopBandFrequency();
        double ap = specs.getPassBandRipple();
        double as = specs.getStopBandAttenuation();

        double wp = 2 * Math.PI * fp;
        double ws = 2 * Math.PI * fs;
        final int n = type.getMinOrderNeeded(ws, wp, ap, as);
        LowPassPrototype lp = new LowPassPrototype(type.buildLowPassPrototype(n, ap, as));
        lp._tf = lpTohp(lp._tf.getNumerator(), lp._tf.getDenominator(), fp);
        lp._tf.normalize();
        return new AnalogFilter(n, lp._tf);
    }

    public static AnalogFilter newBandPass(BandPassSpecs specs, ApproximationType type) {
        double fp1 = specs.getLowerPassBandFrequency();
        double fp2 = specs.getUpperPassBandFrequency();
        double fs1 = specs.getLowerStopBandFrequency();
        double fs2 = specs.getUpperStopBandFrequency();
        double ap = specs.getPassBandRipple();
        double as1 = specs.getLowerStopBandAttenuation();
        double as2 = specs.getUpperStopBandAttenuation();

        double wp1 = 2 * Math.PI * fp1;
        double wp2 = 2 * Math.PI * fp2;
        double ws1 = 2 * Math.PI * fs1;
        double ws2 = 2 * Math.PI * fs2;

        double w0 = Math.sqrt(wp1 * wp2);
        double Q = w0 / (wp2 - wp1);

        double whs1 = ws1 / w0;
        double whs2 = ws2 / w0;

        double omega1 = Q * Math.abs((whs1 * whs1 - 1) / whs1);
        double omega2 = Q * Math.abs((whs2 * whs2 - 1) / whs2);

        final int n1 = type.getMinOrderNeeded(1, omega1, ap, as1);
        final int n2 = type.getMinOrderNeeded(1, omega2, ap, as2);

        int n = 0;
        double as = 0.0;
        if (n1 > n2) {
            n = n1;
            as = as1;
        } else {
            n = n2;
            as = as2;
        }
        TransferFunction tf = zpkToTF(type.buildLowPassPrototype(n, ap, as));

        double bw = fp2 - fp1;
        double f0 = w0 / (2 * Math.PI);
        tf = lpTobp(tf.getNumerator(), tf.getDenominator(), f0, bw);
        tf.normalize();
        return new AnalogFilter(2 * n, tf);
    }

    public static TransferFunction lpTobp(Polynomial num, Polynomial den, double w0, double bw) {
        Polynomial s = new Polynomial(bw, 0.0);
        Polynomial s2w02 = new Polynomial(1.0, 0, w0 * w0);

        RationalFunction bp = new RationalFunction(num, den);
        bp.substituteInPlace(new RationalFunction(s2w02, s));

        return new TransferFunction(bp);
    }

    public static TransferFunction lpTobp(double[] num, double[] den, double w0, double bw) {
        return lpTobp(new Polynomial(num), new Polynomial(den), w0, bw);
    }

    public static AnalogFilter newBandStop(BandStopSpecs specs, ApproximationType type) {
        double fp1 = specs.getLowerPassBandFrequency();
        double fp2 = specs.getUpperPassBandFrequency();
        double fs1 = specs.getLowerStopBandFrequency();
        double fs2 = specs.getUpperStopBandFrequency();
        double ap = specs.getPassBandRipple();
        double as = specs.getStopBandAttenuation();

        double[] wp = new double[2];
        // maximize the pass band
        // https://github.com/scipy/scipy/blob/master/scipy/signal/filter_design.py
        wp[0] = goldenSectionMinimizer(AnalogFilter::bandStopObjMinimize, fp1, fs1 - 1e-12,
                1e-5, 500, specs, type, 0);

        wp[1] = goldenSectionMinimizer(AnalogFilter::bandStopObjMinimize, fs2 + 1e-12, fp2,
                1e-5, 500, specs, type, 1);

        double w1 = (specs.getLowerStopBandFrequency() * (wp[0] - wp[1]) /
                (specs.getLowerStopBandFrequency() * specs.getLowerStopBandFrequency() - wp[0] * wp[1]));
        double w2 = (specs.getUpperStopBandFrequency() * (wp[0] - wp[1]) /
                (specs.getUpperStopBandFrequency() * specs.getUpperStopBandFrequency() - wp[0] * wp[1]));

        double ws = Math.min(Math.abs(w1), Math.abs(w2));
        final int n = type.getMinOrderNeeded(1, ws, ap, as);
        TransferFunction tf = zpkToTF(type.buildLowPassPrototype(n, ap, as));

        double bw = wp[1] - wp[0];
        double f0 = Math.sqrt(wp[0] * wp[1]);
        tf = lpTobs(tf.getNumerator(), tf.getDenominator(), f0, bw);
        tf.normalize();
        return new AnalogFilter(2 * n, tf);
    }

    private static double bandStopObjMinimize(double wp, Object... params) {
        BandStopSpecs specs = (BandStopSpecs) params[0];
        ApproximationType type = (ApproximationType) params[1];
        Integer index = (Integer) params[2];

        double[] passb = new double[]{specs.getLowerPassBandFrequency(), specs.getUpperPassBandFrequency()};
        double[] stopb = new double[]{specs.getLowerStopBandFrequency(), specs.getUpperStopBandFrequency()};
        double amax = specs.getPassBandRipple();
        double amin = specs.getStopBandAttenuation();

        passb[index] = wp;

        double w1 = (stopb[0] * (passb[0] - passb[1]) /
                (stopb[0] * stopb[0] - passb[0] * passb[1]));
        double w2 = (stopb[1] * (passb[0] - passb[1]) /
                (stopb[1] * stopb[1] - passb[0] * passb[1]));

        double w0 = Math.min(Math.abs(w1), Math.abs(w2));
        return type.getOrderNeeded(1, w0, amax, amin);
    }

    public static TransferFunction lpTobs(Polynomial num, Polynomial den, double w0, double bw) {
        Polynomial s = new Polynomial(bw, 0.0);
        Polynomial s2w02 = new Polynomial(1.0, 0.0, w0 * w0);

        RationalFunction bp = new RationalFunction(num, den);
        bp.substituteInPlace(new RationalFunction(s, s2w02));

        return new TransferFunction(bp);
    }

    public static TransferFunction lpTobs(double[] num, double[] den, double w0, double bw) {
        return lpTobs(new Polynomial(num), new Polynomial(den), w0, bw);
    }

    public static TransferFunction lpTohp(Polynomial num, Polynomial den, double w0) {
        final int numDegree = num.degree();
        final int denDegree = den.degree();
        final int filterOrder = Math.max(numDegree, denDegree) + 1;

        // Reverse coefficients then scale them by the
        // order of the denominator i.e. pad with zeros
        double[] hpNumerator = new double[filterOrder];
        for (int i = numDegree, j = 0; i >= 0; --i, ++j) {
            hpNumerator[j] = num.getCoefficientAt(i) * Math.pow(w0, j);
        }

        // Reverse coefficients then scale them by the
        // order of the numerator i.e. pad with zeros
        double[] hpDenominator = new double[filterOrder];
        for (int i = denDegree, j = 0; i >= 0; --i, ++j) {
            hpDenominator[j] = den.getCoefficientAt(i) * Math.pow(w0, j);
        }
        return new TransferFunction(hpNumerator, hpDenominator);
    }

    public static TransferFunction lpTolp(Polynomial num, Polynomial den, double wo) {
        TransferFunction tf = new TransferFunction(num, den);
        tf.substituteInPlace(1.0 / wo);
        return tf;
    }

    public static TransferFunction zpkToTF(ZeroPoleGain zpk) {
        return new TransferFunction(zpk.getZeros(), zpk.getPoles(), zpk.getGain());
    }

    public static void main(String[] args) {
        UnivariateFunction gg = x -> x * x;
        BiFunction<Double, Object[], Double> func = (x, params) -> gg.evaluateAt(x);
        System.out.println(brentsMinimizer(func, 1, 2, 1e-5, 500));
    }
}