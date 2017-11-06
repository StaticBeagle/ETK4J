package com.wildbitsfoundry.etk4j.util;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.MathETK;

public final class NumArrays {
	
	private NumArrays() {}
	
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
	
	public static double[] subtract(double[] a, double[] b) {
		if(a.length != b.length) {
			throw new IllegalArgumentException("a and b dimensions must match");
		}
		final int length = a.length;
		double[] result = new double[length];
		for(int i = 0; i < length; ++i) {
			result[i] = a[i] - b[i];
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
	
	public static double[] multiply(double[] a, double d) {
		double[] result = Arrays.copyOf(a, a.length);
		multiplyInPlace(result, d);
		return result;
	}
	
	public static void multiplyInPlace(double[] a, double d) {
		final int length = a.length;
		for (int i = 0; i < length; ++i) {
			a[i] *= d;
		}
	}
	
	public static double[] multiplyElementWise(double[] a, double[] b) {
		if(a.length != b.length) {
			throw new IllegalArgumentException("a and b dimensions must match");
		}
		
		double[] result = Arrays.copyOf(a, a.length);
		multiplyElementWiseInPlace(result, b);
		return result;
	}
	
	public static void multiplyElementWiseInPlace(double[] a, double[] b) {
		if(a.length != b.length) {
			throw new IllegalArgumentException("a and b dimensions must match");
		}
		
		final int length = a.length;
		for (int i = 0; i < length; ++i) {
			a[i] *= b[i];
		}
	}
	
	public static double[] divideElementWise(double[] a, double[] b) {
		if(a.length != b.length) {
			throw new IllegalArgumentException("a and b dimensions must match");
		}
		
		double[] result = Arrays.copyOf(a, a.length);
		divideElementWiseInPlace(result, b);
		return result;
	}
	
	public static void divideElementWiseInPlace(double[] a, double[] b) {
		if(a.length != b.length) {
			throw new IllegalArgumentException("a and b dimensions must match");
		}
		
		final int length = a.length;
		for (int i = 0; i < length; ++i) {
			a[i] /= b[i];
		}
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
	
	public static double normNegInf(double[] a) {
		double min = Math.abs(a[0]);
		for(int i = 1; i < a.length; ++i) {
			double abs = Math.abs(a[i]);
			if(abs < min) {
				min = abs;
			}
		}
		return min;
	}
	
	/***
	 * Euclidean distance between two arrays
	 * @param a - first array
	 * @param b - second array
	 * @return
	 * 	<pre> sqrt((a[i] - b[i])<sup>2</sup>)
	 */
	public static double distance(double[] a, double[] b) {
		double norm = 0.0;
		for(int i = 0; i < a.length; ++i) {
			norm = MathETK.hypot(norm, a[i] - b[i]);
		}
		return norm;
	}
	
	public static double[] kron(final double[] a, final double[] b) {
		int aLength = a.length;
		int bLength = b.length;
		double[] c = new double[aLength * bLength];
		if (a.length == 1) {
			return NumArrays.multiply(b, a[0]);
		}
		if (b.length == 1) {
			return NumArrays.multiply(a, b[0]);
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
	
	public static double[] repeat(double[] a, int n) {
		final int length = a.length;
		double[] result = new double[length * n];
		for(int i = 0; i < n; ++i) {
			System.arraycopy(a, 0, result, i * length, length);
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
	
	public static double sum(double[] a) {
		final int length = a.length;
		double result = 0;
		for(int i = 0; i < length; ++i) {
			result += a[i];
		}
		return result;
	}
	
	public static double[] cummulativeSum(double[] a) {
		final int length = a.length;
		double[] result = new double[length];
		result[0] = a[0];
		for(int i = 1; i < length; ++i) {
			result[i] = result[i - 1] + a[i];
		}
		return result;
	}
	
	public static double rms(double[] a) {
		final int N = a.length;
		return norm(a) / Math.sqrt(N);
	}
	
	public static boolean isAscending(double[] a) {
		final int n = a.length;
		for(int i = 0; i < n - 1; ++i) {
			if(a[i] > a[i + 1]) {
				return false;
			}
		}
		return true;
	}
	
	/**
	 * Helper method to copy a 2 dimensional array
	 * 
	 * @param array
	 *            to copy
	 * @return a newly created array containing the copy of array
	 */
	public static double[] flatten(double[][] a) {
		int rows = a.length;
		int cols = a[0].length;
		double[] result = new double[rows * cols];
		for (int i = 0; i < rows; i++) {
			System.arraycopy(a[i], 0, result, i * cols, cols);
		}
		return result;
	}

	public static double mean(double[] y) {
		final int n = y.length;
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            sum += y[i];
        }
        return sum / n;
	}
	
	public static Double[] box(double[] a) {
		return Arrays.stream(a).boxed().toArray(size -> new Double[size]);
	}
	
	public static double[] unbox(Double[] a) {
		return Arrays.stream(a).mapToDouble(Double::doubleValue).toArray();
	}
}
