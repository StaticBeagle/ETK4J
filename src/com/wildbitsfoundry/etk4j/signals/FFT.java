package com.wildbitsfoundry.etk4j.signals;

import java.util.Random;

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
		this.inverse(real, imag, FFTScaling.STANDARD);
	}
		
	public void inverse(double[] real, double[] imag, FFTScaling scaling) {
		
		this.direct(imag, real);
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

	public static int uniform(int a, int b) {
		if ((b <= a) || ((long) b - a >= Integer.MAX_VALUE)) {
			throw new IllegalArgumentException("invalid range: [" + a + ", " + b + ")");
		}
		return a + uniform(b - a);
	}

	public static int uniform(int n) {
		if (n <= 0)
			throw new IllegalArgumentException("argument must be positive: " + n);
		return random.nextInt(n);
	}

	private static Random random; // pseudo-random number generator
	private static long seed; // pseudo-random number generator seed

	// static initializer
	static {
		// this is how the seed was set in Java 1.4
		seed = System.currentTimeMillis();
		random = new Random(seed);
	}

	public static double uniform() {
		return random.nextDouble();
	}

	public static double uniform(double a, double b) {
		if (!(a < b)) {
			throw new IllegalArgumentException("invalid range: [" + a + ", " + b + ")");
		}
		return a + uniform() * (b - a);
	}

	public static void main(String[] args) {

		// benchMark();
		int n = 32;
		// double[] real = new double[n];
		// double[] imag = new double[n];
		// for(int i = 0; i < n; ++i) {
		// real[i] = uniform(-1.0, 1.0);
		// imag[i] = uniform(-1.0, 1.0);
		// }

		double[] real = { 0.559677760695923, 0.811961988413724, -0.489763299068234, -0.778067179542341,
				-0.706735964715945, 0.789019721173435, -0.700001555401106, 0.700489895029539, -0.718488129867960,
				0.574342509973089, 0.685516617700664, 0.817398984407037, -0.267529687111015, -0.224422049737993,
				-0.606118364173603, -0.757421636418101, -0.554515146424656, -0.440790757929746, -0.409781781881146,
				0.0141205685354671, 0.280704586278388, -0.452294658802835, -0.681799182902633, 0.864458500245397,
				-0.435253596623234, -0.0249062599816210, 0.177507536424897, 0.918235637008785, 0.659079582021293,
				-0.0474811926324961, -0.285224755579650, 0.863443049844534 };

		double[] imag = { -0.627553075609633, 0.546655823856343, -0.946010653090666, 0.243006473934334,
				-0.257873950960501, 0.253985941117167, 0.549619389302165, 0.278368441048426, -0.896404939976776,
				-0.588221017267259, -0.988610592161680, -0.667865653695379, 0.555211953225667, 0.657655783982322,
				-0.248329940445516, 0.462132700569790, 0.937594247545697, 0.863758857076813, 0.381492588056045,
				-0.151998214979134, -0.929646800306250, 0.0416082390322683, -0.215594991256026, -0.976012577178314,
				0.847837107109384, -0.0561434192544805, 0.531367878780680, -0.296828765390418, -0.0170543373258241,
				-0.222993692923477, 0.775753049706212, -0.142807967492327 };

		for (int i = 0; i < n; ++i) {
			System.out.printf("(%.4f, %.4f)%n", real[i], imag[i]);
		}
		FFT fft = new FFT(n);
		fft.direct(real, imag);

		System.out.println("FFT(x)");
		for (int i = 0; i < n; ++i) {
			System.out.printf("%d: (%.4f, %.4f)%n", i, real[i], imag[i]);
		}

		fft.inverse(real, imag);
		System.out.println("IFFT(x)");
		for (int i = 0; i < n; ++i) {
			System.out.printf("%d: (%.4f, %.4f)%n", i, real[i], imag[i]);
		}
	}

	public static void benchMark() {
		int n = 2048;
		double[] real = new double[n];
		double[] imag = new double[n];
		for (int i = 0; i < n; ++i) {
			real[i] = uniform(-1.0, 1.0);
			imag[i] = uniform(-1.0, 1.0);
		}

		long startTime = System.currentTimeMillis();
		FFT fft = new FFT(n);
		for (int i = 0; i < 100000; ++i) {
			// fft(real, imag);
			fft.direct(real, imag);
		}

		long endTime = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.println(totalTime);
	}
}
