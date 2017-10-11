package com.wildbitsfoundry.etk4j.systems.continuoustime.filters;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.systems.continuoustime.TransferFunction;
import com.wildbitsfoundry.etk4j.util.ArrayUtils;

public class Butterworth extends AnalogFilter {
	
	private TransferFunction _hs = null;
	private int _order;
	private double _eps;
	
	/***
	 * Calculate the minimum order required for Low-Pass Butterworth filter
	 * @param wp passband frequency
	 * @param ws stopband frequency
	 * @param Ap passband attenuation
	 * @param As stopband attenuation
	 * @return
	 */
	public static int getMinOrderRequired(double wp, double ws, double ap, double as) {
		double amax = Math.pow(10, ap * 0.1) - 1;
		double amin = Math.pow(10, as * 0.1) - 1;
		
		double L = Math.log10(amin / amax) / (2 * Math.log10(ws / wp));
		
		return (int)Math.ceil(L);
		
	}
	
	private Butterworth(int n, double ap) {
		_eps = Math.sqrt(Math.pow(10, ap * 0.1) - 1);
		
		final double pid = Math.PI / 180.0;
		Complex[] poles = new Complex[n];
		if(n % 2 == 0) {
			int i = 0;
			for(double k : ArrayUtils.linsteps(-n * 0.5 + 1.0, n * 0.5, 1)) {
				double phik = 180.0 * (k / n) - 90.0 / n;
				poles[i++] = new Complex(-Math.cos(phik * pid), Math.sin(phik * pid));
			}
		} else {
			int i = 0;
			for(double k : ArrayUtils.linsteps(-(n - 1) * 0.5, (n - 1) * 0.5, 1)) {
				double phik = 180.0 * (k / n);
				poles[i++] = new Complex(-Math.cos(phik * pid), Math.sin(phik * pid));
			}
		}
		_hs = new TransferFunction(new Complex[0], poles);
		_order = n;
	}
	
	public static Butterworth newLowPass(double wp, double ws, double ap, double as) {
		final int n = getMinOrderRequired(wp, ws, ap, as);
		Butterworth lp = new Butterworth(n, ap);
		double factor = Math.pow(lp._eps, -1.0 / n) * wp;
		lp._hs.scale(1.0 / factor);
		return lp;
	}
	
	public static Butterworth newLowPass(int n, double ap) {
		return new Butterworth(n, ap);
	}
	
	public static Butterworth newHighPass(double wp, double ws, double ap, double as) {
		final int n = getMinOrderRequired(ws, wp, ap, as);
		Butterworth hp = new Butterworth(n, ap);
		double factor = wp / Math.pow(hp._eps, -1.0 / n);
		hp._hs.scale(factor);
		hp._hs = lpTohp(hp._hs.getNumerator(), hp._hs.getDenominator());
		return hp;
	}
	
	
	public static Butterworth newHighPass(int n, double ap) {
		Butterworth hp = newLowPass(n, ap);
		hp._hs = lpTohp(hp._hs.getNumerator(), hp._hs.getDenominator());
		return hp;
	}
	
	public static TransferFunction lpTohp(Polynomial numerator, Polynomial denominator) {
		return lpTohp(numerator.getCoefficients(), denominator.getCoefficients());
	}
	
	public static TransferFunction lpTohp(double[] numerator, double[] denominator) {
		final int numDegree = numerator.length - 1;
		final int denDegree = denominator.length - 1;
		final int filterOrder = numDegree + denDegree + 1;
		
		// Reverse coefficients then scale them by the 
		// order of the denominator i.e. pad with zeros
		double[] hpNumerator = new double[filterOrder];
		for(int i = numDegree, j = 0; i >= 0; --i, ++j) {
			hpNumerator[j] = numerator[i];
		}
		
		// Reverse coefficients then scale them by the 
		// order of the numerator i.e. pad with zeros
		double[] hpDenominator = new double[filterOrder];
		for(int i = denDegree, j = 0; i >= 0; --i, ++j) {
			hpDenominator[j] = denominator[i];
		}
		
		return new TransferFunction(hpNumerator, hpDenominator);
	}
	
	public double getEpsilon() {
		return _eps;
	}
	
	public int getOrder() {
		return _order;
	}
	
//	public TransferFunction getTF() {
//		return _hs;
//	}
	
	public static void main(String[] args) {
		Butterworth lowpass = Butterworth.newLowPass(1 * 2 * Math.PI, 10 * 2 * Math.PI, 0.2, 60);
		
		Butterworth highpass = Butterworth.newHighPass(10 * 2 * Math.PI, 1 * 2 * Math.PI, 0.2, 60);
		
		Polynomial f = new Polynomial(new double[] {1, 1});
		System.out.println(f.pow(2));
		System.out.println(f.pow(3));
		System.out.println(f.substitute(new Polynomial(new double[] {1, 1})));
		System.out.println();
		
		System.out.printf("Low pass: %n%s%n%n", lowpass._hs.toString());
		System.out.printf("High pass: %n%s%n", highpass._hs.toString());
	}
}
