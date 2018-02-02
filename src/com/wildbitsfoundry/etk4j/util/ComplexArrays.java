package com.wildbitsfoundry.etk4j.util;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

public final class ComplexArrays {
	private ComplexArrays() {

	}

	public static Complex[] deepCopy(Complex[] a) {
		final int length = a.length;
		Complex[] result = new Complex[length];
		for (int i = 0; i < length; ++i) {
			result[i] = new Complex(a[i]);
		}
		return result;
	}

	public static Complex[] conv(Complex[] a, Complex[] b) {
		Complex[] result = new Complex[a.length + b.length - 1];
		for (int i = 0; i < result.length; ++i) {
			result[i] = new Complex();
			for (int j = Math.max(0, i + 1 - b.length); j < Math.min(a.length, i + 1); ++j) {
				result[i].addEquals(a[j].multiply(b[i - j]));
			}
		}
		return result;
	}
	
	/***
	 * Element wise 
	 * @param a
	 * @return
	 */
	public static double[] real(Complex[] a) {
		double[] result = new double[a.length];
		for(int i = 0; i < a.length; ++i) {
			result[i] = a[i].real();
		}
		return result;
	}
	
	public static double[] imag(Complex[] a) {
		double[] result = new double[a.length];
		for(int i = 0; i < a.length; ++i) {
			result[i] = a[i].imag();
		}
		return result;
	}
	
	/***
	 * 
	 * @param a
	 * @return
	 */
	public static Complex[] conj(Complex[] a) {
		Complex[] result = new Complex[a.length];
		for(int i = 0; i < a.length; ++i) {
			result[i] = a[i].conj();
		}
		return result;
	}
}