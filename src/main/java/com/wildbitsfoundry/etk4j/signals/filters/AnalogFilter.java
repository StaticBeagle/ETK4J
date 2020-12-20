package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.LowPassSpecs;

import java.util.function.BiFunction;

public class AnalogFilter {

    private int _order;
    private TransferFunction _tf;

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


    /***
     * Calculate the minimum order required for Low-Pass Chebyshev filter
     *
     * @param fp
     *            passband frequency in Hertz
     * @param fs
     *            stopband frequency in Hertz
     * @param ap
     *            passband attenuation in dB
     * @param as
     *            stopband attenuation in dB
     * @return
     */
    public static int getMinOrderNeeded(double fp, double fs, double ap, double as, ApproximationType type) {
        return type.getMinOrderNeeded(fp, fs, ap, as);
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
        return new AnalogFilter(n, tf);
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
		wp[0] = goldenSectionMinimum(AnalogFilter::bandStopObjMinimize, fp1, fs1 - 1e-12,
				1e-5, 500, specs, type, 0);

        wp[1] = goldenSectionMinimum(AnalogFilter::bandStopObjMinimize, fs2 + 1e-12, fp2,
                1e-5, 500, specs, type, 1);

        double w1 = (specs.getLowerStopBandFrequency() * (wp[0] - wp[1]) /
                (specs.getLowerStopBandFrequency() * specs.getLowerStopBandFrequency() - wp[0] * wp[1]));
        double w2 = (specs.getUpperStopBandFrequency() * (wp[0] - wp[1]) /
                (specs.getUpperStopBandFrequency() * specs.getUpperStopBandFrequency() - wp[0] * wp[1]));

        double ws = Math.min(Math.abs(w1), Math.abs(w2));
        final int n = type.getMinOrderNeeded(1, ws, ap, as);
        TransferFunction tf = zpkToTF(type.buildLowPassPrototype(n, ap, as));

        double bw = fp2 - fp1;
        double f0 = Math.sqrt(fp1 * fp2);
        tf = lpTobs(tf.getNumerator(), tf.getDenominator(), f0, bw);
        tf.normalize();
        return new AnalogFilter(n, tf);
    }

    // TODO move this somewhere where it fits better
    // Maybe improve it and roll out our own implementation
    // https://www.mathworks.com/matlabcentral/fileexchange/25919-golden-section-method-algorithm
	private static double goldenSectionMinimum(BiFunction<Double, Object[], Double> func,
                                               double a, double b, double tol, int maxIter, Object... params) {
		double gold = (Math.sqrt(5.0) - 1.0) / 2.0;

		double x1 = a + (1 - gold) * (b - a);
		double x2 = a + gold * (b - a);

		double fx1 = func.apply(x1, params);
		double fx2 = func.apply(x2, params);

		int k = 1;
		while ((Math.abs(b - a) > tol) && (k < maxIter)) {
			k = k + 1;
			if (fx1 < fx2) {
				b = x2;
				x2 = x1;
				x1 = a + (1 - gold) * (b - a);
			} else {
				a = x1;
				x1 = x2;
				x2 = a + gold * (b - a);
			}
			fx1 = func.apply(x1, params);
			fx2 = func.apply(x2, params);
		}
		if (fx1 < fx2)
			return x1;
		return x2;
	}

    private static double bandStopObjMinimize(double wp, Object ... params) {
    	BandStopSpecs specs = (BandStopSpecs) params[0];
    	ApproximationType type = (ApproximationType) params[1];
    	Integer index = (Integer) params[2];

    	double[] passb = new double[] {specs.getLowerPassBandFrequency(), specs.getUpperPassBandFrequency()};
    	double[] stopb = new double[] {specs.getLowerStopBandFrequency(), specs.getUpperStopBandFrequency()};
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

    public int getOrder() {
        return _order;
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
        return new TransferFunction(zpk.Zeros, zpk.Poles, zpk.Gain);
    }
}