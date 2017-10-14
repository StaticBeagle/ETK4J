package com.wildbitsfoundry.etk4j.util;

import java.util.Arrays;

public final class ArrayUtils {
	
//	 /**
//	  * Updates MACHEPS for the executing machine.
//	  */
//	 public static void updateMacheps() {
//	  MACHEPS = 1;
//	  do
//	   MACHEPS /= 2;
//	  while (1 + MACHEPS / 2 != 1);
//	 }
	private ArrayUtils() {}
	
	/***
	 * Creates n linearly spaced samples between x0 and x1
	 * @param x0 starting point
	 * @param x1 end point
	 * @param n number of samples
	 * @return Array of n equally spaced samples between [x0, x1]
	 */
	public static double[] linspace(double x0, double x1, int n) {
		double[] result = new double[n];
		int i = 0;
		double delta = (x1 - x0) / --n;
		while(i < n) {
			result[i] = x0 + i * delta;
			++i;
		}
		result[i] = x1;
		return result;
	}
	
	/***
	 * Creates n linearly spaced samples between x0 and x1
	 * @param x0 starting point
	 * @param x1 end point
	 * @param step step size
	 * @return Array of n equally spaced samples between [x0, x1]
	 */
	public static double[] linsteps(double x0, double x1, double step) {
		final int n = (int) Math.ceil((x1 - x0) / step);
		double[] result = new double[n + 1];
		int i = 0;
		while(i < n) {
			result[i] = x0 + i * step;
			++i;
		}
		result[i] = x1;
		return result;
	}
	
	/***
	 * Creates n logarithmically spaced samples between decades x0 and x1
	 * @param x0 starting decade
	 * @param x1 end decade
	 * @param n number of samples
	 * @return Array of n logarithmically spaced samples between [x0, x1]
	 */
	public static double[] logspace(int x0, int x1, int n) {
		double[] result = new double[n];
		int i = 0;
		double delta = (double)(x1 - x0) / --n;
		while(i < n) {
			result[i] = Math.pow(10.0, x0 + i * delta);
			++i;
		}
		result[i] = Math.pow(10.0, x1);
		return result;
	}
	
	public static double[] add(double[] a, double[] b) {
		if(a.length != b.length) {
			throw new IllegalArgumentException("a and b dimensions must match");
		}
		final int length = a.length;
		double[] result = new double[length];
		for(int i = 0; i < length; ++i) {
			result[i] = a[i] + b[i];
		}
		return result;
	}
	
	public static double[] conv(double[] a, double[] b) {
		double[] result = new double[a.length + b.length - 1];

		for (int i = 0; i < result.length; ++i) {
			result[i] = 0.0;
			for (int j = Math.max(0, i + 1 - b.length); j < Math.min(a.length, i + 1); ++j) {
				result[i] += a[j] * b[i - j];
			}
		}
		return result;
	}
	
	public static double[] conv(double[] a, double d) {
		double[] result = new double[a.length];

		for (int i = 0; i < result.length; ++i) {
			result[i] = a[i] * d;
		}
		return result;
	}
	
	public static double normFast(double[] a) {
		double norm = 0.0;
		for(int i = 0; i < a.length; ++i) {
			norm += a[i] * a[i];
		}
		return Math.sqrt(norm);
	}
	
	public static double max(double[] a) {
		double max = a[0];
		for(int i = 1; i < a.length; ++i) {
			if(a[i] > max) {
				max = a[i];
			}
		}
		return max;
	}
	
	public static double norm1(double[] a) {
		double norm = 0.0;
		for(int i = 0; i < a.length; ++i) {
			norm += Math.abs(a[i]);
		}
		return norm;
	}
	
	public static double norm(double[] a) {
		double max = normInf(a);
		if(max == 0.0) {
			return 0.0;
		}
		double tmp = Math.pow(1.0 / max, 2);
		double norm = 0.0;
		for(int i = 0; i < a.length; ++i) {
			norm += a[i] * a[i] * tmp;
		}
		return Math.sqrt(norm) * max;
	}
	
	public static double normInf(double[] a) {
		double max = Math.abs(a[0]);
		for(int i = 1; i < a.length; ++i) {
			double abs = Math.abs(a[i]);
			if(abs > max) {
				max = abs;
			}
		}
		return max;
	}
	
	public static double[] kron(final double[] a, final double[] b) {
		int aLength = a.length;
		int bLength = b.length;
		double[] c = new double[aLength * bLength];
		if (a.length == 1) {
			return ArrayUtils.conv(b, a[0]);
		}
		if (b.length == 1) {
			return ArrayUtils.conv(a, b[0]);
		}
		for (int i = 0; i < aLength; i++) {
			for (int j = 0; j < bLength; j++) {
				c[j + i * (bLength)] = a[i] * b[j];
			}
		}
		return c;
	}
	
	public static double[] concat(double[] a, double[] b) {
		int aLength = a.length;
		int bLength = b.length;
		double[] result = new double[aLength + bLength];
		System.arraycopy(a, 0, result, 0, aLength);
		System.arraycopy(b, 0, result, aLength, bLength);
		return result;
	}

	public static double[] concatAll(double[] a, double[]... rest) {
		int totalLength = a.length;
		for (double[] array : rest) {
			totalLength += array.length;
		}
		double[] result = Arrays.copyOf(a, totalLength);
		int offset = a.length;
		for (double[] array : rest) {
			System.arraycopy(array, 0, result, offset, array.length);
			offset += array.length;
		}
		return result;
	}
	
	public static double[] reverse(double[] a) {
		final int length = a.length;
		double[] result = new double[a.length];
		for(int i = length - 1, j = 0; i >= 0; --i, ++j) {
			result[j] = a[i];
		}
		return result;
	}
}
