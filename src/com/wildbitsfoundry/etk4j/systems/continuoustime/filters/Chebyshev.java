package com.wildbitsfoundry.etk4j.systems.continuoustime.filters;

import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.systems.continuoustime.TransferFunction;
import com.wildbitsfoundry.etk4j.util.NumArrays;

public class Chebyshev extends AnalogFilter {

	private TransferFunction _tf = null;
	private double _eps;

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
	public static int getMinOrderRequired(double fp, double fs, double ap, double as) {
		double wp = 2 * Math.PI * fp;
		double ws = 2 * Math.PI * fs;
		double amax = Math.pow(10, ap * 0.1) - 1;
		double amin = Math.pow(10, as * 0.1) - 1;

		double L = MathETK.acosh(Math.sqrt(amin / amax)) / MathETK.acosh(ws / wp);

		return (int) Math.ceil(L);
	}

	private Chebyshev(int n, double ap) {
		_eps = Math.sqrt(Math.pow(10, ap * 0.1) - 1);

		double a = 1.0 / n * MathETK.asinh(1 / _eps);
		double sinha = Math.sinh(a);
		double cosha = Math.cosh(a);

		final double pid = Math.PI / 180.0;
		Complex[] poles = new Complex[n];
		if (n % 2 == 0) {
			int i = 0;
			for (double k : NumArrays.linsteps(-n * 0.5 + 1.0, n * 0.5, 1)) {
				double phik = 180.0 * (k / n) - 90.0 / n;
				poles[i++] = new Complex(-sinha * Math.cos(phik * pid), cosha * Math.sin(phik * pid));
			}
		} else {
			int i = 0;
			for (double k : NumArrays.linsteps(-(n - 1) * 0.5, (n - 1) * 0.5, 1)) {
				double phik = 180.0 * (k / n);
				poles[i++] = new Complex(-sinha * Math.cos(phik * pid), cosha * Math.sin(phik * pid));
			}
		}

		double N = 1;
		for (int k = 0; k < (int) Math.ceil(n / 2.0); ++k) {
			if (poles[k].imag() != 0.0) {
				N *= poles[k].real() * poles[k].real() + poles[k].imag() * poles[k].imag();
			} else {
				N *= -poles[k].real();
			}
		}
		_tf = new TransferFunction(N, poles);
		_order = n;
	}

	public static Chebyshev newLowPass(int n, double ap) {
		return new Chebyshev(n, ap);
	}

	public static Chebyshev newLowPass(double fp, double fs, double ap, double as) {
		double wp = 2 * Math.PI * fp;
		double ws = 2 * Math.PI * fs;
		final int n = getMinOrderRequired(wp, ws, ap, as);
		Chebyshev lp = new Chebyshev(n, ap);
		double factor = Math.pow(lp._eps, -1.0 / n) * wp;
		lp._tf.substituteInPlace(1.0 / factor);
		return lp;
	}

	public static Chebyshev newHighPass(int n, double ap) {
		Chebyshev hp = newLowPass(n, ap);
		hp._tf = lpTohp(hp._tf.getNumerator(), hp._tf.getDenominator());
		return hp;
	}

	public static Chebyshev newHighPass(double fp, double fs, double ap, double as) {
		double wp = 2 * Math.PI * fp;
		double ws = 2 * Math.PI * fs;
		final int n = getMinOrderRequired(ws, wp, ap, as);
		Chebyshev hp = new Chebyshev(n, ap);
		double factor = wp / Math.pow(hp._eps, -1.0 / n);
		hp._tf.substituteInPlace(factor);
		hp._tf = lpTohp(hp._tf.getNumerator(), hp._tf.getDenominator());
		return hp;
	}

	public static Chebyshev newBandPass(double fp1, double fp2, double fs1, double fs2, double ap, double as1,
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
		Chebyshev bp = new Chebyshev(n, ap);
		double bw = Q / Math.pow(bp._eps, -1.0 / n) / w0;
		bp._tf = lpTobp(bp._tf.getNumerator(), bp._tf.getDenominator(), w0, bw);
		bp._order <<= 1;
		return bp;
	}

	public static Chebyshev newBandStop(double fp1, double fp2, double fs1, double fs2, double amax, double amin) {
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
		Chebyshev bp = new Chebyshev(n, amax);
		double bw = Q * Math.pow(bp._eps, -1.0 / n) / w0;
		bp._tf = lpTobs(bp._tf.getNumerator(), bp._tf.getDenominator(), w0, bw);
		bp._order <<= 1;
		return bp;
	}

	public static void main(String[] args) {
		Chebyshev lowpass = newLowPass(1, 2, 1, 20);

		Chebyshev highpass = newHighPass(10, 1, 0.2, 60);

		Chebyshev bandpass = newBandPass(190, 210, 180, 220, 0.2, 20, 20);

		Chebyshev bandstop = newBandStop(3.6e3, 9.1e3, 5.45e3, 5.90e3, 1.5, 38);

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
