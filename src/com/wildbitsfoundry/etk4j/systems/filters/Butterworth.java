package com.wildbitsfoundry.etk4j.systems.filters;

import static com.wildbitsfoundry.etk4j.systems.filters.AnalogFilter.lpTobp;
import static com.wildbitsfoundry.etk4j.systems.filters.AnalogFilter.lpTobs;
import static com.wildbitsfoundry.etk4j.systems.filters.AnalogFilter.lpTohp;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.wildbitsfoundry.etk4j.constants.ETKConstants;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.systems.TransferFunction;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.BandPassSpecs;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.LowPassSpecs;
import com.wildbitsfoundry.etk4j.util.NumArrays;
import static com.wildbitsfoundry.etk4j.math.specialfunctions.Elliptic.compEllipInt1;

public class Butterworth {

	private TransferFunction _tf = null;
	private double _eps;

	/***
	 * Calculate the minimum order required for Low-Pass Butterworth filter
	 * 
	 * @param wp
	 *            passband frequency
	 * @param ws
	 *            stopband frequency
	 * @param Ap
	 *            passband attenuation
	 * @param As
	 *            stopband attenuation
	 * @return
	 */
	public static int getMinOrderRequired(double fp, double fs, double ap, double as) {
		double wp = 2 * Math.PI * fp;
		double ws = 2 * Math.PI * fs;
		double amax = Math.pow(10, ap * 0.1) - 1;
		double amin = Math.pow(10, as * 0.1) - 1;

		double L = Math.log10(amin / amax) / (2 * Math.log10(ws / wp));

		return (int) Math.ceil(L);

	}

	private Butterworth(int n, double ap) {
		_eps = Math.sqrt(Math.pow(10, ap * 0.1) - 1);

		final double pid = Math.PI / 180.0;
		Complex[] poles = new Complex[n];
		if (n % 2 == 0) {
			int i = 0;
			for (double k : NumArrays.linsteps(-n * 0.5 + 1.0, 1, n * 0.5)) {
				double phik = 180.0 * (k / n) - 90.0 / n;
				poles[i++] = new Complex(-Math.cos(phik * pid), Math.sin(phik * pid));
			}
		} else {
			int i = 0;
			for (double k : NumArrays.linsteps(-(n - 1) * 0.5, 1, (n - 1) * 0.5)) {
				double phik = 180.0 * (k / n);
				poles[i++] = new Complex(-Math.cos(phik * pid), Math.sin(phik * pid));
			}
		}
		_tf = new TransferFunction(new Complex[0], poles);
	}

	public static Butterworth newLowPass(int n, double ap) {
		return new Butterworth(n, ap);
	}

	public static Butterworth newLowPass(double fp, double fs, double ap, double as) {
		double wp = 2 * Math.PI * fp;
		double ws = 2 * Math.PI * fs;
		final int n = getMinOrderRequired(wp, ws, ap, as);
		Butterworth lp = new Butterworth(n, ap);
		double factor = Math.pow(lp._eps, -1.0 / n) * wp;
		lp._tf.substituteInPlace(1.0 / factor);
		return lp;
	}

	public static Butterworth newHighPass(int n, double ap) {
		Butterworth hp = newLowPass(n, ap);
		hp._tf = lpTohp(hp._tf.getNumerator(), hp._tf.getDenominator(), 0.0);
		return hp;
	}

	public static Butterworth newHighPass(double fp, double fs, double ap, double as) {
		double wp = 2 * Math.PI * fp;
		double ws = 2 * Math.PI * fs;
		final int n = getMinOrderRequired(ws, wp, ap, as);
		Butterworth hp = new Butterworth(n, ap);
		double factor = wp / Math.pow(hp._eps, -1.0 / n);
		hp._tf.substituteInPlace(factor);
		hp._tf = lpTohp(hp._tf.getNumerator(), hp._tf.getDenominator(), 0.0);
		return hp;
	}

	public static Butterworth newBandPass(double fp1, double fp2, double fs1, double fs2, double ap, double as1,
			double as2) {
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

		final int n1 = getMinOrderRequired(1, omega1, ap, as1);
		final int n2 = getMinOrderRequired(1, omega2, ap, as2);

		final int n = Math.max(n1, n2);
		Butterworth bp = new Butterworth(n, ap);
		double bw = Q / Math.pow(bp._eps, -1.0 / n) / w0;
		bp._tf = lpTobp(bp._tf.getNumerator(), bp._tf.getDenominator(), w0, bw);
		return bp;
	}

	public static Butterworth newBandStop(double fp1, double fp2, double fs1, double fs2, double amax, double amin) {
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

		final int n1 = getMinOrderRequired(1, omegas1, amax, amin);
		final int n2 = getMinOrderRequired(1, omegas2, amax, amin);

		final int n = Math.max(n1, n2);
		Butterworth bp = new Butterworth(n, amax);
		double bw = Q * Math.pow(bp._eps, -1.0 / n) / w0;
		bp._tf = lpTobs(bp._tf.getNumerator(), bp._tf.getDenominator(), w0, bw);
		return bp;
	}

	public static void main(String[] args) {
		testButter(); 	// working 
		testCheby1();	// working
	}

	public static void testButter() {
		LowPassSpecs lpSpecs = new LowPassSpecs();
		lpSpecs.PassBandAttenuation = -20 * Math.log10(1.0 / Math.sqrt(2.0));
		lpSpecs.StopBandAttenuation = 70;
		lpSpecs.PassBandFrequency = 1.0 / (2 * Math.PI);
		lpSpecs.StopBandFrequency = 10.0 / (2 * Math.PI);
		AnalogFilter lp = AnalogFilter.newLowPass(lpSpecs, ApproximationType.BUTTERWORTH);
		// Move these two to unit tests they are ready
		System.out.println(Arrays.toString(lp._tf.getNumerator().getCoefficients()));
		System.out.println(Arrays.toString(lp._tf.getDenominator().getCoefficients()));

		HighPassSpecs hpSpecs = new HighPassSpecs();
		hpSpecs.PassBandAttenuation = -20 * Math.log10(1.0 / Math.sqrt(2.0));
		hpSpecs.StopBandAttenuation = 70;
		hpSpecs.PassBandFrequency = 1.0 / (2 * Math.PI);
		hpSpecs.StopBandFrequency = 0.1 / (2 * Math.PI);
		AnalogFilter hp = AnalogFilter.newHighPass(hpSpecs, ApproximationType.BUTTERWORTH);
		// Move these two to unit tests they are ready
		System.out.println(Arrays.toString(hp._tf.getNumerator().getCoefficients()));
		System.out.println(Arrays.toString(hp._tf.getDenominator().getCoefficients()));

		BandPassSpecs bpSpecs = new BandPassSpecs();
		bpSpecs.LowerPassBandFrequency = 190;
		bpSpecs.UpperPassBandFrequency = 210;
		bpSpecs.LowerStopBandFrequency = 180;
		bpSpecs.UpperStopBandFrequency = 220;
		bpSpecs.PassBandAttenuation = -20 * Math.log10(1.0 / Math.sqrt(2.0));
		bpSpecs.LowerStopBandAttenuation = 20;
		bpSpecs.UpperStopBandAttenuation = 20;
		AnalogFilter bp = AnalogFilter.newBandPass(bpSpecs, ApproximationType.BUTTERWORTH);
		// Move these two to unit tests they are ready
		System.out.println(Arrays.toString(bp._tf.getNumerator().getCoefficients()));
		System.out.println(Arrays.toString(bp._tf.getDenominator().getCoefficients()));

		BandStopSpecs bsSpecs = new BandStopSpecs();
		bsSpecs.LowerPassBandFrequency = 3.6e3;
		bsSpecs.UpperPassBandFrequency = 9.1e3;
		bsSpecs.LowerStopBandFrequency = 5.45e3;
		bsSpecs.UpperStopBandFrequency = 5.90e3;
		bsSpecs.PassBandAttenuation = -20 * Math.log10(1.0 / Math.sqrt(2.0));
		bsSpecs.StopBandAttenuation = 38;
		AnalogFilter bs = AnalogFilter.newBandStop(bsSpecs, ApproximationType.BUTTERWORTH);
		// Move these two to unit tests they are ready
		System.out.println(Arrays.toString(bs._tf.getNumerator().getCoefficients()));
		System.out.println(Arrays.toString(bs._tf.getDenominator().getCoefficients()));

		System.out.printf("Low pass: %n%s%n%n", lp._tf.toString());
		System.out.printf("High pass: %n%s%n%n", hp._tf.toString());
		System.out.printf("Band pass: %n%s%n%n", bp._tf.toString());
		System.out.printf("Band stop: %n%s%n", bs._tf.toString());
	}
	
	public static void testCheby1() {
		LowPassSpecs lpSpecs = new LowPassSpecs();
		lpSpecs.PassBandAttenuation = 0.5;
		lpSpecs.StopBandAttenuation = 70;
		lpSpecs.PassBandFrequency = 1.0 / (2 * Math.PI);
		lpSpecs.StopBandFrequency = 10.0 / (2 * Math.PI);
		AnalogFilter lp = AnalogFilter.newLowPass(lpSpecs, ApproximationType.CHEBYSHEV);
		// Move these two to unit tests they are ready
		System.out.println(Arrays.toString(lp._tf.getNumerator().getCoefficients()));
		System.out.println(Arrays.toString(lp._tf.getDenominator().getCoefficients()));

		HighPassSpecs hpSpecs = new HighPassSpecs();
		hpSpecs.PassBandAttenuation = 0.5;
		hpSpecs.StopBandAttenuation = 70;
		hpSpecs.PassBandFrequency = 1.0 / (2 * Math.PI);
		hpSpecs.StopBandFrequency = 0.1 / (2 * Math.PI);
		AnalogFilter hp = AnalogFilter.newHighPass(hpSpecs, ApproximationType.CHEBYSHEV);
		// Move these two to unit tests they are ready
		System.out.println(Arrays.toString(hp._tf.getNumerator().getCoefficients()));
		System.out.println(Arrays.toString(hp._tf.getDenominator().getCoefficients()));

		BandPassSpecs bpSpecs = new BandPassSpecs();
		bpSpecs.LowerPassBandFrequency = 190;
		bpSpecs.UpperPassBandFrequency = 210;
		bpSpecs.LowerStopBandFrequency = 180;
		bpSpecs.UpperStopBandFrequency = 220;
		bpSpecs.PassBandAttenuation = 0.5;
		bpSpecs.LowerStopBandAttenuation = 20;
		bpSpecs.UpperStopBandAttenuation = 20;
		AnalogFilter bp = AnalogFilter.newBandPass(bpSpecs, ApproximationType.CHEBYSHEV);
		// Move these two to unit tests they are ready
		System.out.println(Arrays.toString(bp._tf.getNumerator().getCoefficients()));
		System.out.println(Arrays.toString(bp._tf.getDenominator().getCoefficients()));

		BandStopSpecs bsSpecs = new BandStopSpecs();
		bsSpecs.LowerPassBandFrequency = 3.6e3;
		bsSpecs.UpperPassBandFrequency = 9.1e3;
		bsSpecs.LowerStopBandFrequency = 5.45e3;
		bsSpecs.UpperStopBandFrequency = 5.90e3;
		bsSpecs.PassBandAttenuation = 0.5;
		bsSpecs.StopBandAttenuation = 38;
		AnalogFilter bs = AnalogFilter.newBandStop(bsSpecs, ApproximationType.CHEBYSHEV);
		// Move these two to unit tests they are ready
		System.out.println(Arrays.toString(bs._tf.getNumerator().getCoefficients()));
		System.out.println(Arrays.toString(bs._tf.getDenominator().getCoefficients()));

		System.out.printf("Low pass: %n%s%n%n", lp._tf.toString());
		System.out.printf("High pass: %n%s%n%n", hp._tf.toString());
		System.out.printf("Band pass: %n%s%n%n", bp._tf.toString());
		System.out.printf("Band stop: %n%s%n", bs._tf.toString());
	}
}
