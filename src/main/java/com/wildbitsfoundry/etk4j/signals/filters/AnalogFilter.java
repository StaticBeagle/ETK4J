package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.LowPassSpecs;

public class AnalogFilter {

	private int _order;
	private TransferFunction _tf;


	public AnalogFilter(LowPassSpecs specs, LowPassPrototype lowPassPrototype) {
		double fp = specs.getPassBandFrequency();
		double fs = specs.getStopBandFrequency();
		double ap = specs.getPassBandRipple();
		double as = specs.getStopBandAttenuation();

		double wp = 2 * Math.PI * fp;
		double ws = 2 * Math.PI * fs;
		final int n = lowPassPrototype.getMinOrderNeeded(wp, ws, ap, as);
		TransferFunction tf = new TransferFunction(lowPassPrototype.buildLowPassPrototype(n, ap, as));
		double w0 = lowPassPrototype.getScalingFrequency(wp, ws);

		this._order = n;
		this._tf = lpTolp(tf.getNumerator(), tf.getDenominator(), w0);

	}

	public AnalogFilter(HighPassSpecs specs, LowPassPrototype lowPassPrototype) {
		double fp = specs.getPassBandFrequency();
		double fs = specs.getStopBandFrequency();
		double ap = specs.getPassBandRipple();
		double as = specs.getStopBandAttenuation();

		double wp = 2 * Math.PI * fp;
		double ws = 2 * Math.PI * fs;
		final int n = lowPassPrototype.getMinOrderNeeded(ws, wp, ap, as);
		TransferFunction tf = new TransferFunction(lowPassPrototype.buildLowPassPrototype(n, ap, as));
		double w0 = lowPassPrototype.getScalingFrequency(wp, ws);

		this._order = n;
		this._tf = lpTohp(tf.getNumerator(), tf.getDenominator(), w0);
	}

	AnalogFilter(BandPassSpecs specs, LowPassPrototype lowPassPrototype) {
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

		final int n1 = lowPassPrototype.getMinOrderNeeded(1, omega1, ap, as1);
		final int n2 = lowPassPrototype.getMinOrderNeeded(1, omega2, ap, as2);

		int n = 0;
		double as = 0.0;
		if (n1 > n2) {
			n = n1;
			as = as1;
		} else {
			n = n2;
			as = as2;
		}
		TransferFunction tf = new TransferFunction(lowPassPrototype.buildLowPassPrototype(n, ap, as));

		double bw = Q / w0;

		this._order = n;
		this._tf = lpTobp(tf.getNumerator(), tf.getDenominator(), w0, bw);
	}

	AnalogFilter(BandStopSpecs specs, LowPassPrototype lowPassPrototype) {
		double fp1 = specs.getLowerPassBandFrequency();
		double fp2 = specs.getUpperPassBandFrequency();
		double fs1 = specs.getLowerStopBandFrequency();
		double fs2 = specs.getUpperStopBandFrequency();
		double amax = specs.getPassBandRipple();
		double amin = specs.getStopBandAttenuation();

		double wp1 = 2 * Math.PI * fp1;
		double wp2 = 2 * Math.PI * fp2;
		double ws1 = 2 * Math.PI * fs1;
		double ws2 = 2 * Math.PI * fs2;
		double w0 = Math.sqrt(wp1 * wp2);
		double Q = w0 / (wp2 - wp1);

		double whs1 = ws1 / w0;
		double whs2 = ws2 / w0;

		double omegas1 = 1 / (Q * Math.abs((whs1 * whs1 - 1) / whs1));
		double omegas2 = 1 / (Q * Math.abs((whs2 * whs2 - 1) / whs2));

		final int n1 = lowPassPrototype.getMinOrderNeeded(1, omegas1, amax, amin);
		final int n2 = lowPassPrototype.getMinOrderNeeded(1, omegas2, amax, amin);

		int n = n1 > n2 ? n1 : n2;

		TransferFunction tf = new TransferFunction(lowPassPrototype.buildLowPassPrototype(n, amax, amin));

		double bw = w0 / Q;
		tf = lpTobs(tf.getNumerator(), tf.getDenominator(), w0, bw);

		this._order = n;
		this._tf = tf;
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

	public static TransferFunction lpTobp(Polynomial num, Polynomial den, double w0, double bw) {
		Polynomial s = new Polynomial(1 / bw, 0.0);
		Polynomial s2w02 = new Polynomial(1.0, 0, w0 * w0);

		RationalFunction bp = new RationalFunction(num, den);
		bp.substituteInPlace(new RationalFunction(s2w02, s));

		return new TransferFunction(bp);
	}

	public static TransferFunction lpTobs(Polynomial num, Polynomial den, double w0, double bw) {
		Polynomial s = new Polynomial(bw, 0.0);
		Polynomial s2w02 = new Polynomial(1.0, 0.0, w0 * w0);

		RationalFunction bp = new RationalFunction(num, den);
		bp.substituteInPlace(new RationalFunction(s, s2w02));

		return new TransferFunction(bp);
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
}
