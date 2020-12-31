package com.wildbitsfoundry.etk4j.signals.fft;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.util.NumArrays;

public class FFT {

	public enum FFTScaling {
		STANDARD,
		UNITARY,
		NONE
	}

	private int _n;
	private int _m;

	private double[] _cos;
	private double[] _sin;

	public FFT(int n) {
		_n = n;
		_m = (int) (Math.log(n) / Math.log(2.0));
		
		if(_n != (1 << _m)) {
			throw new IllegalArgumentException("n must be a power of 2");
		}

		_cos = new double[n / 2];
		_sin = new double[n / 2];

		double t = -2 * Math.PI / n;
		for (int i = 0; i < n / 2; ++i) {
			_cos[i] = Math.cos(i * t);
			_sin[i] = Math.sin(i * t);
		}
	}

	/***************************************************************
	 * fft.c Douglas L. Jones University of Illinois at Urbana-Champaign January
	 * 19, 1992 http://cnx.rice.edu/content/m12016/latest/
	 * 
	 * fft: in-place radix-2 DIT DFT of a complex input
	 * 
	 * input: n: length of FFT: must be a power of two m: n = 2**m input/output
	 * x: double array of length n with real part of data y: double array of
	 * length n with imag part of data
	 * 
	 * Permission to copy and use this program is granted as long as this header
	 * is included.
	 ****************************************************************/
	public void direct(double[] real, double[] imag) {
		if(real.length != imag.length) {
			throw new IllegalArgumentException("Length mismatch between real and imag");
		}
		
		int i, j, k, n1, n2, a;
		double c, s, t1, t2;

		// Bit-reverse
		j = 0;
		n2 = _n / 2;
		for (i = 1; i < _n - 1; i++) {
			n1 = n2;
			while (j >= n1) {
				j = j - n1;
				n1 = n1 / 2;
			}
			j = j + n1;

			if (i < j) {
				t1 = real[i];
				real[i] = real[j];
				real[j] = t1;
				t1 = imag[i];
				imag[i] = imag[j];
				imag[j] = t1;
			}
		}

		// FFT
		n1 = 0;
		n2 = 1;

		for (i = 0; i < _m; i++) {
			n1 = n2;
			n2 = n2 + n2;
			a = 0;

			for (j = 0; j < n1; j++) {
				c = _cos[a];
				s = _sin[a];
				a += 1 << (_m - i - 1);

				for (k = j; k < _n; k = k + n2) {
					t1 = c * real[k + n1] - s * imag[k + n1];
					t2 = s * real[k + n1] + c * imag[k + n1];
					real[k + n1] = real[k] - t1;
					imag[k + n1] = imag[k] - t2;
					real[k] = real[k] + t1;
					imag[k] = imag[k] + t2;
				}
			}
		}
	}

	public void inverse(double[] real, double[] imag) {
		this.direct(imag, real);
		this.inverse(real, imag, FFTScaling.STANDARD);
	}

	public void inverse(double[] real, double[] imag, FFTScaling scaling) {


		switch(scaling) {
		case STANDARD:
			double factor = 1.0 / _n;
			NumArrays.multiplyInPlace(real, factor);
			NumArrays.multiplyInPlace(imag, factor);
			break;
		case UNITARY:
			factor = 1.0 / Math.sqrt(_n);
			NumArrays.multiplyInPlace(real, factor);
			NumArrays.multiplyInPlace(imag, factor);
		case NONE:
			break;
		default:
			throw new IllegalStateException("FFT inverse DFTScaling missing");
		}
	}
}
