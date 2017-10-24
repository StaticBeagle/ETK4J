package com.wildbitsfoundry.etk4j.signals;

import java.util.Random;

public class FFT {

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
	public static void fft(double[] x, double[] y) {
		int i, j, k, n1, n2, a;
		double c, s, t1, t2;

		int n = x.length;
		int m = (int) (Math.log(n) / Math.log(2.0));

		// need to do lazy evaluation and cache this
		// tables to optimize this
		double[] cos = new double[n / 2];
		double[] sin = new double[n / 2];

		for (int ii = 0; ii < n / 2; ii++) {
			cos[ii] = Math.cos(-2 * Math.PI * ii / n);
			sin[ii] = Math.sin(-2 * Math.PI * ii / n);
		}

		// Bit-reverse
		j = 0;
		n2 = n / 2;
		for (i = 1; i < n - 1; i++) {
			n1 = n2;
			while (j >= n1) {
				j = j - n1;
				n1 = n1 / 2;
			}
			j = j + n1;

			if (i < j) {
				t1 = x[i];
				x[i] = x[j];
				x[j] = t1;
				t1 = y[i];
				y[i] = y[j];
				y[j] = t1;
			}
		}

		// FFT
		n1 = 0;
		n2 = 1;

		for (i = 0; i < m; i++) {
			n1 = n2;
			n2 = n2 + n2;
			a = 0;

			for (j = 0; j < n1; j++) {
				c = cos[a];
				s = sin[a];
				a += 1 << (m - i - 1);

				for (k = j; k < n; k = k + n2) {
					t1 = c * x[k + n1] - s * y[k + n1];
					t2 = s * x[k + n1] + c * y[k + n1];
					x[k + n1] = x[k] - t1;
					y[k + n1] = y[k] - t2;
					x[k] = x[k] + t1;
					y[k] = y[k] + t2;
				}
			}
		}
	}

	private static final double[] W_SUB_N_R = { 0x1.0p0, -0x1.0p0, 0x1.1a62633145c07p-54, 0x1.6a09e667f3bcdp-1,
			0x1.d906bcf328d46p-1, 0x1.f6297cff75cbp-1, 0x1.fd88da3d12526p-1, 0x1.ff621e3796d7ep-1, 0x1.ffd886084cd0dp-1,
			0x1.fff62169b92dbp-1, 0x1.fffd8858e8a92p-1, 0x1.ffff621621d02p-1, 0x1.ffffd88586ee6p-1,
			0x1.fffff62161a34p-1, 0x1.fffffd8858675p-1, 0x1.ffffff621619cp-1, 0x1.ffffffd885867p-1,
			0x1.fffffff62161ap-1, 0x1.fffffffd88586p-1, 0x1.ffffffff62162p-1, 0x1.ffffffffd8858p-1,
			0x1.fffffffff6216p-1, 0x1.fffffffffd886p-1, 0x1.ffffffffff621p-1, 0x1.ffffffffffd88p-1,
			0x1.fffffffffff62p-1, 0x1.fffffffffffd9p-1, 0x1.ffffffffffff6p-1, 0x1.ffffffffffffep-1,
			0x1.fffffffffffffp-1, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0,
			0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0,
			0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0, 0x1.0p0,
			0x1.0p0 };

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

		int n = 32;
		// double[] real = new double[n];
		// double[] imag = new double[n];
		// for(int i = 0; i < n; ++i) {
		// real[i] = uniform(-1.0, 1.0);
		// imag[i] = uniform(-1.0, 1.0);
		// }

		double[] real = new double[] { 0.559677760695923, 0.811961988413724, -0.489763299068234, -0.778067179542341,
				-0.706735964715945, 0.789019721173435, -0.700001555401106, 0.700489895029539, -0.718488129867960,
				0.574342509973089, 0.685516617700664, 0.817398984407037, -0.267529687111015, -0.224422049737993,
				-0.606118364173603, -0.757421636418101, -0.554515146424656, -0.440790757929746, -0.409781781881146,
				0.0141205685354671, 0.280704586278388, -0.452294658802835, -0.681799182902633, 0.864458500245397,
				-0.435253596623234, -0.0249062599816210, 0.177507536424897, 0.918235637008785, 0.659079582021293,
				-0.0474811926324961, -0.285224755579650, 0.863443049844534 };

		double[] imag = new double[] { -0.627553075609633, 0.546655823856343, -0.946010653090666, 0.243006473934334,
				-0.257873950960501, 0.253985941117167, 0.549619389302165, 0.278368441048426, -0.896404939976776,
				-0.588221017267259, -0.988610592161680, -0.667865653695379, 0.555211953225667, 0.657655783982322,
				-0.248329940445516, 0.462132700569790, 0.937594247545697, 0.863758857076813, 0.381492588056045,
				-0.151998214979134, -0.929646800306250, 0.0416082390322683, -0.215594991256026, -0.976012577178314,
				0.847837107109384, -0.0561434192544805, 0.531367878780680, -0.296828765390418, -0.0170543373258241,
				-0.222993692923477, 0.775753049706212, -0.142807967492327 };

		for (int i = 0; i < n; ++i) {
			System.out.printf("(%.4f, %.4f)%n", real[i], imag[i]);
		}

		fft(real, imag);

		System.out.println("FFT(x)");
		for (int i = 0; i < n; ++i) {
			System.out.printf("%d: (%.4f, %.4f)%n", i, real[i], imag[i]);
		}
	}
}
