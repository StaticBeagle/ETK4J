package com.wildbitsfoundry.etk4j.systems.continuoustime.filters;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.systems.continuoustime.TransferFunction;
import com.wildbitsfoundry.etk4j.util.NumArrays;

public class Butterworth extends AnalogFilter {

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
			for (double k : NumArrays.linsteps(-n * 0.5 + 1.0, n * 0.5, 1)) {
				double phik = 180.0 * (k / n) - 90.0 / n;
				poles[i++] = new Complex(-Math.cos(phik * pid), Math.sin(phik * pid));
			}
		} else {
			int i = 0;
			for (double k : NumArrays.linsteps(-(n - 1) * 0.5, (n - 1) * 0.5, 1)) {
				double phik = 180.0 * (k / n);
				poles[i++] = new Complex(-Math.cos(phik * pid), Math.sin(phik * pid));
			}
		}
		_tf = new TransferFunction(new Complex[0], poles);
		_order = n;
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
		hp._tf = lpTohp(hp._tf.getNumerator(), hp._tf.getDenominator());
		return hp;
	}

	public static Butterworth newHighPass(double fp, double fs, double ap, double as) {
		double wp = 2 * Math.PI * fp;
		double ws = 2 * Math.PI * fs;
		final int n = getMinOrderRequired(ws, wp, ap, as);
		Butterworth hp = new Butterworth(n, ap);
		double factor = wp / Math.pow(hp._eps, -1.0 / n);
		hp._tf.substituteInPlace(factor);
		hp._tf = lpTohp(hp._tf.getNumerator(), hp._tf.getDenominator());
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
		bp._order <<= 1;
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
		bp._order <<= 1;
		return bp;
	}

	public static void main(String[] args) {
		Butterworth lowpass = newLowPass(1, 10, 0.2, 60);

		Butterworth highpass = newHighPass(10, 1, 0.2, 60);

		Butterworth bandpass = newBandPass(190, 210, 180, 220, 0.2, 20, 20);

		Butterworth bandstop = newBandStop(3.6e3, 9.1e3, 5.45e3, 5.90e3, 1.5, 38);

		Polynomial f = new Polynomial(new double[] { 1, 1 });
		System.out.println(f.pow(2));
		System.out.println(f.pow(3));
		System.out.println(f.substitute(new Polynomial(new double[] { 1, 1 })));
		System.out.println();

		System.out.printf("Low pass: %n%s%n%n", lowpass._tf.toString());
		System.out.printf("High pass: %n%s%n%n", highpass._tf.toString());
		System.out.printf("Band pass: %n%s%n%n", bandpass._tf.toString());
		System.out.printf("Band stop: %n%s%n", bandstop._tf.toString());
	}
}
