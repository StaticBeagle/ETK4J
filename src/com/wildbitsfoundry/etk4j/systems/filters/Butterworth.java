package com.wildbitsfoundry.etk4j.systems.filters;

import static com.wildbitsfoundry.etk4j.systems.filters.AnalogFilter.lpTobp;
import static com.wildbitsfoundry.etk4j.systems.filters.AnalogFilter.lpTobs;
import static com.wildbitsfoundry.etk4j.systems.filters.AnalogFilter.lpTohp;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.systems.TransferFunction;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.BandPassSpecs;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.LowPassSpecs;
import com.wildbitsfoundry.etk4j.util.NumArrays;

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
	
	// K(k) = RF(0, 1 - k^2), 1);
	public static double CompleteEllipticalIntegral(double k, double tol) {
		return RF(0.0, 1 - k * k, 1.0, tol);
	}
		
	// Carlson's elliptical integral of the first kind
	public static double RF(double x, double y, double z, double tol) {
		// check bounds
		
		double dx, dy, dz;
		double lambda;
		double n = 1.0 / 3.0;
		double mean;
		double tmp;
		do {
			lambda = Math.sqrt(x * y) + Math.sqrt(y * z) + Math.sqrt(z * x);
			x = 0.25 * (x + lambda);
			y = 0.25 * (y + lambda);
			z = 0.25 * (z + lambda);
			mean = (x + y + z) * n;
			tmp = 1 / mean;
			dx = (mean - x) * tmp;
			dy = (mean - y) * tmp;
			dz = (mean - z) * tmp;
		}while(Math.max(Math.max(Math.abs(dx), Math.abs(dy)), Math.abs(dz)) > tol);
		double e2 = dx * dy - dz * dz;
		double e3 = dx * dy * dz;
		double c1 = 1.0 / 24.0;
		double c2 = 0.1;
		double c3 = 3.0 / 44.0;
		double c4 = 1.0 / 14.0;
		
		double result = 1.0 + (c1 * e2 - c2 - c3 * e3) * e2 + c4 * e3;
		return result / Math.sqrt(mean);
	}

	public static void main(String[] args) {
		LowPassSpecs lpSpecs = new LowPassSpecs();
		lpSpecs.PassBandAttenuation = 1;
		lpSpecs.StopBandAttenuation = 25;
		lpSpecs.PassBandFrequency = 0.6 / (2.0 * Math.PI);
		lpSpecs.StopBandFrequency = 1.0 / (2.0 * Math.PI);
		AnalogFilter lp = AnalogFilter.newLowPass(lpSpecs, ApproximationType.INVERSE_CHEBYSHEV);

		HighPassSpecs hpSpecs = new HighPassSpecs();
		hpSpecs.PassBandAttenuation = 0.2;
		hpSpecs.StopBandAttenuation = 60;
		hpSpecs.PassBandFrequency = 10;
		hpSpecs.StopBandFrequency = 1;
		AnalogFilter hp = AnalogFilter.newHighPass(hpSpecs, ApproximationType.BUTTERWORTH);

		BandPassSpecs bpSpecs = new BandPassSpecs();
		bpSpecs.LowerPassBandFrequency = 190;
		bpSpecs.UpperPassBandFrequency = 210;
		bpSpecs.LowerStopBandFrequency = 180;
		bpSpecs.UpperStopBandFrequency = 220;
		bpSpecs.PassBandAttenuation = 0.2;
		bpSpecs.LowerStopBandAttenuation = 20;
		bpSpecs.UpperStopBandAttenuation = 20;
		AnalogFilter bp = AnalogFilter.newBandPass(bpSpecs, ApproximationType.BUTTERWORTH);

		BandStopSpecs bsSpecs = new BandStopSpecs();
		bsSpecs.LowerPassBandFrequency = 3.6e3;
		bsSpecs.UpperPassBandFrequency = 9.1e3;
		bsSpecs.LowerStopBandFrequency = 5.45e3;
		bsSpecs.UpperStopBandFrequency = 5.90e3;
		bsSpecs.PassBandAttenuation = 1.5;
		bsSpecs.StopBandAttenuation = 38;
		AnalogFilter bs = AnalogFilter.newBandStop(bsSpecs, ApproximationType.INVERSE_CHEBYSHEV);

		System.out.printf("Low pass: %n%s%n%n", lp._tf.toString());
		System.out.printf("High pass: %n%s%n%n", hp._tf.toString());
		System.out.printf("Band pass: %n%s%n%n", bp._tf.toString());
		System.out.printf("Band stop: %n%s%n", bs._tf.toString());

		long start = System.currentTimeMillis();
		double total = 0.0;
		double current = 0.0;
		for (int i = 0; i < 10; ++i) {
			for (int ii = 0; ii < 1000000; ii++) {
				//current = Math.sqrt(5 * 5 + 5 * 5);
				current = MathETK.hypot(5d, 5d);
			}
			total += current;
		}
		long elapsed = System.currentTimeMillis() - start;
		System.out.println("elapsed time = " + elapsed + "ms");
		System.out.println();
		
		
		double CEI = CompleteEllipticalIntegral(0.5, 0.08);
		System.out.println(CEI);
		
		
		double ap = 1;
		double as = 34;
		double wp = 1;
		double ws = 2;
		
		double k = wp / ws;
		double kp = Math.pow(1 - k * k, 0.25);
		double eps = Math.sqrt(Math.pow(10.0, 0.1 * ap) - 1);
		double kn =  eps / Math.sqrt(Math.pow(10.0, 0.1 * as) - 1);
		
		double ne = CompleteEllipticalIntegral(k, 0.08) * CompleteEllipticalIntegral(Math.sqrt(1 - kn * kn), 0.08);
		ne /= CompleteEllipticalIntegral(Math.sqrt(1 - k * k), 0.08) * CompleteEllipticalIntegral(kn, 0.08);
		


		double u = 0.5 * (1 - kp) / (1 + kp);
		double q = u + 2 * Math.pow(u, 5) + 15 * Math.pow(u, 9) + 150 * Math.pow(u, 13);
		
		double D = (Math.pow(10.0, 0.1 * as) - 1) / (Math.pow(10.0, 0.1 * ap) - 1);
		
		double nexact = (Math.log10(16.0 * D) / Math.log10(1.0 / q));
		int n = (int) Math.ceil(nexact);
		System.out.println(n);
		double V = 0.5 / n * Math.log((Math.pow(10.0, ap * 0.05) + 1) / (Math.pow(10.0, ap * 0.05) - 1));
		
		double p0num = 0.0;
		double p0den = 0.0;
		for(int i = 0; i < 5;){
			p0num += Math.pow(-1.0, i) * Math.pow(q, i * (i + 1.0)) * Math.sinh((2.0 * i + 1) * V);
			++i;
			p0den += Math.pow(-1.0, i) * Math.pow(q, i * i) * Math.cosh(2 * i * V);
		}
		p0num *= Math.pow(q, 0.25);
		p0den += 0.5;
		double p0 = Math.abs(p0num / p0den);
		
		double p02 = p0 * p0;
		double W = Math.sqrt((1 + p02 / k) * (1 + k * p02));
		
		int r = n - (n % 2) >> 1;
		double[] a = new double[r];
		double[] b = new double[r];
		double[] c = new double[r];
		for(int i = 0; i < r; ++i) {
			double xinum = 0.0;
			double xiden = 0.0;
			double mu = n % 2 == 0 ? i + 0.5 : i + 1.0;
			double mup = mu * Math.PI / n;
			for(int m = 0; m < 5;){
				xinum += Math.pow(-1.0, m) * Math.pow(q, m * (m + 1.0)) * Math.sin((2.0 * m + 1) * mup);
				++m;
				xiden += Math.pow(-1.0, m) * Math.pow(q, m * m) * Math.cos(2 * m * mup);
			}
			xinum *= 2.0 * Math.pow(q, 0.25);
			xiden = 1.0 + 2.0 * xiden;
			double xi = xinum / xiden;
			double Xi2 = xi  * xi;
			double yi = Math.sqrt((1.0 - Xi2 / k) * (1 - k * Xi2));
			double tmp = 1.0 / (1 + p02 * Xi2);
			a[i] = 1.0 / Xi2;
			b[i] = (2 * p0 * yi) * tmp;
			c[i] = (p02 * yi * yi + Xi2 * W * W) * tmp * tmp;
		}
		
		double H0 = 1.0;
		for(int i = 0; i < r; ++i) {
			H0 *= c[i] / a[i];
		}
		H0 *= n % 2 == 0 ? Math.pow(10.0, -ap / 20.0) : p0;
		Polynomial num = new Polynomial(1, 0, a[0]);
		Polynomial den = new Polynomial(1, b[0], c[0]);
		for(int i = 1; i < r; ++i) {
			num.multiplyEquals(1, 0, a[i]);
			den.multiplyEquals(1, b[i], c[i]);
		}
		num.multiplyEquals(H0);
		if(n % 2 != 0) {
			den.multiplyEquals(1, p0);
		}
		System.out.println(p0);
		
		TransferFunction tf = new TransferFunction(num, den);
		System.out.println(tf);
		
		System.out.println(Arrays.toString(a));
		System.out.println(Arrays.toString(b));
		System.out.println(Arrays.toString(c));
		
		System.out.println(Arrays.toString(tf.getNumerator().getCoefficients()));
		System.out.println(Arrays.toString(tf.getDenominator().getCoefficients()));
	}
}
