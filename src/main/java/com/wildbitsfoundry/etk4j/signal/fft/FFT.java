package com.wildbitsfoundry.etk4j.signal.fft;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

/**
 * The {@code FFT} class provides and implementation of the Fast Fourier Transform.
 * @see <a href="https://en.wikipedia.org/wiki/Fast_Fourier_transform">Fast Fourier Transform.</a>
 */
public class FFT {

	private int _n;
	private int _m;

	private double[] _cos;
	private double[] _sin;

	/**
	 * Construcsts and instance of the {@code FFT} class.
	 * @param n The length of the {@code FFT}. Must be a power of 2.
	 */
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
	/**
	 * Fast Fourier Transform in place.The real and imaginary parts after performing the {@code FFT}, are stored in the
	 * input arrays {@code real} and {@code imag} respectively.
	 * @param real The real part of the data.
	 * @param imag The imaginary part of the data.
	 */
	public void direct(double[] real, double[] imag) {
		if(real.length != imag.length) {
			throw new IllegalArgumentException("Length mismatch between real and imag");
		}
		if(real.length != _n) {
			throw new IllegalArgumentException(String.format("The lengths of the arrays must be equal to n = %d.", _n));
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

	/**
	 * Inverse Fast Fourier Transform in place.
	 * @param real The real part of the data.
	 * @param imag The imaginary part of the data.
	 */
	public void inverse(double[] real, double[] imag) {
		this.direct(imag, real);
		double factor = 1.0 / _n;
		multiplyInPlace(real, imag, factor);
	}
	
    private static void multiplyInPlace(double[] a, double[] b, double d) {
        final int length = a.length;
        for (int i = 0; i < length; ++i) {
            a[i] *= d;
            b[i] *= d;
        }
    }

	/**
	 * Fast Fourier Transform in place.
	 * @param data The input data.
	 */
	public void direct(Complex[] data) {
		if(data.length != _n) {
			throw new IllegalArgumentException(String.format("The lengths of the arrays must be equal to n = %d.", _n));
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
				// swap data[i] and data[j]
				Complex swap = data[i];
				data[i] = data[j];
				data[j] = swap;
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
					t1 = c * data[k + n1].real() - s * data[k + n1].imag();
					t2 = s * data[k + n1].real() + c * data[k + n1].imag();
					data[k + n1].addEquals(-data[k + n1].real() + data[k].real() - t1, -data[k + n1].imag() + data[k].imag() - t2);
					data[k].addEquals(t1, t2);
				}
			}
		}
	}

	/**
	 * Inverse Fast Fourier Transform in place.
	 * @param data The input data.
	 */
	public void inverse(Complex[] data) {
		double t1;
		// finding the conjugate of the complex number
		// without calling conj() and creating a new Complex
		for(int i = 0; i < _n; ++i) {
			t1 = data[i].imag();
			data[i].addEquals(0.0, -t1 - t1);
		}
		this.direct(data);
		double factor = 1.0 / _n;
		for(int i = 0; i < _n; ++i) {
			t1 = data[i].imag();
			data[i].addEquals(0.0, -t1 - t1);
			data[i].multiplyEquals(factor);
		}
	}
}
