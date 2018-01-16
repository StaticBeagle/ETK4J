package com.wildbitsfoundry.etk4j.systems.filters;

import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;
import com.wildbitsfoundry.etk4j.systems.TransferFunction;
import com.wildbitsfoundry.etk4j.systems.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.BandPassSpecs;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.LowPassSpecs;

public class AnalogFilter {
	
	int _order;
	TransferFunction _tf;
	
	static class LowPassPrototype {
		private TransferFunction _tf;
		
		public LowPassPrototype(TransferFunction tf) {
			_tf = tf;
		}
		
		public LowPassPrototype(ZeroPoleGain zpk) {
			_tf = zpkToTF(zpk);
		}
	}
	
	protected AnalogFilter(int order, TransferFunction tf) {
		_order = order;
		_tf = tf;
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
	
//	protected static AnalogFilter newLowPass(int n, double ap, ApproximationType type) {
//		TransferFunction lp = zpkToTF(type.buildLowPassPrototype(n, ap, 0.0));
//		return new AnalogFilter(n, lp);
//	}
	
	public static AnalogFilter newLowPass(LowPassSpecs specs, ApproximationType type) {
		double fp = specs.PassBandFrequency;
		double fs = specs.StopBandFrequency;
		double ap = specs.PassBandAttenuation;
		double as = specs.StopBandAttenuation;
		
		double wp = 2 * Math.PI * fp;
		double ws = 2 * Math.PI * fs;
		final int n = type.getMinOrderNeeded(wp, ws, ap, as);
		double lpa = type.getAttenuation(ap, as);
		LowPassPrototype lp = new LowPassPrototype(type.buildLowPassPrototype(n, lpa, as));
		double w0 = type.getScalingFrequency(wp,ws);
		lp._tf = lpTolp(lp._tf.getNumerator(), lp._tf.getDenominator(), w0);
		return new AnalogFilter(n, lp._tf);
	}
	
//	public static AnalogFilter newHighPass(int n, double ap, ApproximationType type) {
//		LowPassPrototype lp = type.buildLowPassPrototype(n, ap, 0.0);
//		lp._tf = lpTohp(lp._tf.getNumerator(), lp._tf.getDenominator(), 0.0);
//		return new AnalogFilter(n, lp._tf);
//	}
	
	public static AnalogFilter newHighPass(HighPassSpecs specs, ApproximationType type) {
		double fp = specs.PassBandFrequency;
		double fs = specs.StopBandFrequency;
		double ap = specs.PassBandAttenuation;
		double as = specs.StopBandAttenuation;
		
		double wp = 2 * Math.PI * fp;
		double ws = 2 * Math.PI * fs;
		final int n = type.getMinOrderNeeded(ws, wp, ap, as);
		double hpa = type.getAttenuation(ap, as);
		LowPassPrototype lp = new LowPassPrototype(type.buildLowPassPrototype(n, hpa, as));
		double factor = type.getScalingFrequency(wp, ws);
		lp._tf = lpTohp(lp._tf.getNumerator(), lp._tf.getDenominator(), factor);
		return new AnalogFilter(n, lp._tf);
	}
	
	public static AnalogFilter newBandPass(BandPassSpecs specs, ApproximationType type) {
		double fp1 = specs.LowerPassBandFrequency;
		double fp2 = specs.UpperPassBandFrequency;
		double fs1 = specs.LowerStopBandFrequency;
		double fs2 = specs.UpperStopBandFrequency;
		double ap = specs.PassBandAttenuation;
		double as1 = specs.LowerStopBandAttenuation;
		double as2 = specs.UpperStopBandAttenuation;
		
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
		if(n1 > n2){
			n = n1;
			as = as1;
		} else {
			n = n2;
			as = as2;
		}
		double bpa = type.getAttenuation(ap, as);
		TransferFunction tf = zpkToTF(type.buildLowPassPrototype(n, bpa, as));

		double bw = Q / w0;
		tf = lpTobp(tf.getNumerator(), tf.getDenominator(), w0, bw);
		return new AnalogFilter(n, tf);
	}
	
	public static AnalogFilter newBandStop(BandStopSpecs specs, ApproximationType type) {
		double fp1 = specs.LowerPassBandFrequency;
		double fp2 = specs.UpperPassBandFrequency;
		double fs1 = specs.LowerStopBandFrequency;
		double fs2 = specs.UpperStopBandFrequency;
		double amax = specs.PassBandAttenuation;
		double amin = specs.StopBandAttenuation;
		
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
		
		final int n1 = type.getMinOrderNeeded(1, omegas1, amax, amin);
		final int n2 = type.getMinOrderNeeded(1, omegas2, amax, amin);
		
		int n = n1 > n2 ? n1 : n2;
		
		double bsa = type.getAttenuation(amax, amin);
		TransferFunction tf = zpkToTF(type.buildLowPassPrototype(n, bsa, amin));

		double bw = w0 / Q;
		tf = lpTobs(tf.getNumerator(), tf.getDenominator(), w0, bw);
		return new AnalogFilter(n, tf);
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
		return  new  TransferFunction(hpNumerator, hpDenominator);
	}
	
	public static TransferFunction lpTolp(Polynomial num, Polynomial den, double wo) {
		TransferFunction tf = new TransferFunction(num, den);
		tf.substituteInPlace(1.0 / wo);
		return tf;
	}
	
	public static TransferFunction zpkToTF(ZeroPoleGain zpk) {
		TransferFunction tf = new TransferFunction(zpk.Zeros, zpk.Poles);
		return tf.multiply(zpk.Gain);
	}
}
