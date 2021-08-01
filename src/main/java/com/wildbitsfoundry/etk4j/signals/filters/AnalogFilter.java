package com.wildbitsfoundry.etk4j.signals.filters;

import static com.wildbitsfoundry.etk4j.math.optimize.minimizers.GoldenSection.goldenSectionMinimizer;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.LowPassSpecs;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;

public class AnalogFilter {

    private static int getRelativeDegree(Complex[] zeros, Complex[] poles) {
        int degree = poles.length - zeros.length;
        if(degree < 0) {
            // throw
        }
        return degree;
    }

    public static TransferFunction lpTolp(ZeroPoleGain zpk, double w0) {
        Complex[] zeros = zpk.getZeros();
        Complex[] poles = zpk.getPoles();
        double k = zpk.getGain();
        int degree = getRelativeDegree(zeros, poles);
        ComplexArrays.multiplyInPlace(zeros, w0);
        ComplexArrays.multiplyInPlace(poles, w0);
        k *= Math.pow(w0, degree);
        return new TransferFunction(zeros, poles, k);
    }

    public static TransferFunction lpTobp(ZeroPoleGain zpk, double w0, double bw) {
        Complex[] zeros = zpk.getZeros();
        Complex[] poles = zpk.getPoles();
        double k = zpk.getGain();
        Complex[] zlp = ComplexArrays.multiply(zeros, bw * 0.5);
        Complex[] plp = ComplexArrays.multiply(poles, bw * 0.5);

        int degree = getRelativeDegree(zeros, poles);

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
        return new TransferFunction(zbp, pbp, k);
    }

    public static TransferFunction lpTohp(ZeroPoleGain zpk, double w0) {
        Complex[] zeros = zpk.getZeros();
        Complex[] poles = zpk.getPoles();
        double k = zpk.getGain();
        int degree = getRelativeDegree(zeros, poles);

        Complex[] zhp = ComplexArrays.divide(w0, zeros);
        Complex[] php = ComplexArrays.divide(w0, poles);

        zhp = ComplexArrays.concat(zhp, ComplexArrays.zeros(degree));

        zeros = Arrays.stream(zeros).map(Complex::uminus).toArray(Complex[]::new);
        poles = Arrays.stream(poles).map(Complex::uminus).toArray(Complex[]::new);
        k *= ComplexArrays.product(zeros).divide(ComplexArrays.product(poles)).real();

        RationalFunction rf = new RationalFunction(zhp, php, k);
        return new TransferFunction(zhp, php, k);
    }

    public static TransferFunction lpTobs(ZeroPoleGain zpk, double w0, double bw) {
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

        return new TransferFunction(zbs, pbs, k);
    }

    protected static FilterOrderResults.OrderAndCutoffFrequency lowPassFilterOrder(LowPassSpecs specs,
                                                                                   FilterOrderCalculationStrategy strategy) {
        double wp = specs.getPassBandFrequency();
        double ws = specs.getStopBandFrequency();
        double rp = specs.getPassBandRipple();
        double rs = specs.getStopBandAttenuation();

        double nat = ws / wp;

        double gStop = Math.pow(10, 0.1 * rs);
        double gPass = Math.pow(10, 0.1 * rp);
        int n = strategy.calculateMinOrder(nat, gPass, gStop);
        double wn = strategy.calculateLowPassWn(n, specs);

        return new FilterOrderResults.OrderAndCutoffFrequency(n, wn);
    }

    protected static FilterOrderResults.OrderAndCutoffFrequency highPassFilterOrder(HighPassSpecs specs, FilterOrderCalculationStrategy strategy) {
        double wp = specs.getPassBandFrequency();
        double ws = specs.getStopBandFrequency();
        double rp = specs.getPassBandRipple();
        double rs = specs.getStopBandAttenuation();

        double nat = wp / ws;

        double gStop = Math.pow(10, 0.1 * rs);
        double gPass = Math.pow(10, 0.1 * rp);
        int n = strategy.calculateMinOrder(nat, gPass, gStop);
        double wn = strategy.calculateHighPassWn(n, specs);

        return new FilterOrderResults.OrderAndCutoffFrequency(n, wn);
    }

    protected static FilterOrderResults.OrderAndCutoffFrequencies bandPassFilterOrder(FilterSpecs.BandPassSpecs specs,
                                                                                FilterOrderCalculationStrategy strategy) {
        double wp1 = specs.getLowerPassBandFrequency();
        double wp2 = specs.getUpperPassBandFrequency();
        double ws1 = specs.getLowerStopBandFrequency();
        double ws2 = specs.getUpperStopBandFrequency();
        double rp = specs.getPassBandRipple();
        double rs = specs.getStopBandAttenuation();

        double w1 = (ws1 * ws1 - wp1 * wp2) / (ws1 * (wp1 - wp2));
        double w2 = (ws2 * ws2 - wp1 * wp2) / (ws2 * (wp1 - wp2));

        double nat = Math.min(Math.abs(w1), Math.abs(w2));

        double gStop = Math.pow(10, 0.1 * rs);
        double gPass = Math.pow(10, 0.1 * rp);

        int n = strategy.calculateMinOrder(nat, gPass, gStop);
        double[] wn = strategy.calculateBandPassWn(n, specs);
        return new FilterOrderResults.OrderAndCutoffFrequencies(n, wn[0], wn[1]);
    }

    protected static FilterOrderResults.OrderAndCutoffFrequencies bandStopFilterOrder(FilterSpecs.BandStopSpecs specs,
                                                                                FilterOrderCalculationStrategy strategy) {
        double wp1 = specs.getLowerPassBandFrequency();
        double wp2 = specs.getUpperPassBandFrequency();
        double ws1 = specs.getLowerStopBandFrequency();
        double ws2 = specs.getUpperStopBandFrequency();
        double rp = specs.getPassBandRipple();
        double rs = specs.getStopBandAttenuation();

        double[] wp = new double[2];
        // maximize the pass band
        // https://github.com/scipy/scipy/blob/master/scipy/signal/filter_design.py
        wp[0] = goldenSectionMinimizer(AnalogFilter::bandStopObjMinimize, wp1, ws1 - 1e-12,
                1e-5, 500, specs, strategy, 0);

        wp[1] = goldenSectionMinimizer(AnalogFilter::bandStopObjMinimize, ws2 + 1e-12, wp2,
                1e-5, 500, specs, strategy, 1);

        wp1 = wp[0];
        wp2 = wp[1];
        double w1 = ws1 * (wp[0] - wp[1]) / (ws1 * ws1 - wp[0] * wp[1]);
        double w2 = ws2 * (wp[0] - wp[1]) / (ws2 * ws2 - wp[0] * wp[1]);

        double nat = Math.min(Math.abs(w1), Math.abs(w2));

        double gStop = Math.pow(10, 0.1 * rs);
        double gPass = Math.pow(10, 0.1 * rp);
        int n = strategy.calculateMinOrder(nat, gPass, gStop);

        BandStopSpecs specsCopy = new BandStopSpecs(specs);
        specsCopy.setLowerPassBandFrequency(wp1);
        specsCopy.setUpperPassBandFrequency(wp2);
        double[] wn = strategy.calculateBandStopWn(n, specsCopy);

        Arrays.sort(wn);
        return new FilterOrderResults.OrderAndCutoffFrequencies(n, wn[0], wn[1]);
    }

    private static double bandStopObjMinimize(double wp, Object... params) {
        BandStopSpecs specs = (BandStopSpecs) params[0];
        FilterOrderCalculationStrategy type = (FilterOrderCalculationStrategy) params[1];
        Integer index = (Integer) params[2];

        double[] passb = new double[]{specs.getLowerPassBandFrequency(), specs.getUpperPassBandFrequency()};
        double[] stopb = new double[]{specs.getLowerStopBandFrequency(), specs.getUpperStopBandFrequency()};
        double rp = specs.getPassBandRipple();
        double rs = specs.getStopBandAttenuation();

        passb[index] = wp;

        double w1 = (stopb[0] * (passb[0] - passb[1]) /
                (stopb[0] * stopb[0] - passb[0] * passb[1]));
        double w2 = (stopb[1] * (passb[0] - passb[1]) /
                (stopb[1] * stopb[1] - passb[0] * passb[1]));

        double gStop = Math.pow(10, 0.1 * rs);
        double gPass = Math.pow(10, 0.1 * rp);
        double nat = Math.min(Math.abs(w1), Math.abs(w2));
        return type.calculateExactOrder(nat, gPass, gStop);
    }
}
