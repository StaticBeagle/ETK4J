package com.wildbitsfoundry.etk4j.util;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.IntStream;

import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkXYDimensions;

public final class ComplexArrays {
	private ComplexArrays() {

	}

	public static Complex[] deepCopy(Complex[] a) {
		final int length = a.length;
		Complex[] result = new Complex[length];
		for (int i = 0; i < length; ++i) {
			result[i] = Complex.newComplex(a[i]);
		}
		return result;
	}

	public static Complex[] convolution(Complex[] a, Complex[] b) {
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

	public static Complex[] fromReal(double[] a) {
		Complex[] result = new Complex[a.length];
		for(int i = 0; i < a.length; ++i) {
			result[i] = Complex.fromReal(a[i]);
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

	public static Complex[] fromImaginary(double[] a) {
		Complex[] result = new Complex[a.length];
		for(int i = 0; i < a.length; ++i) {
			result[i] = Complex.fromImaginary(a[i]);
		}
		return result;
	}

	public static double[] abs(Complex[] a) {
		if(a.length == 0) {
			return new double[0];
		}
		if(a == null) {
			return null;
		}
		double[] result = new double[a.length];
		for(int i = 0; i < a.length; ++i) {
			result[i] = a[i].abs();
		}
		return result;
	}
	
	public static Complex[] zip(double[] real, double[] imag) {
		checkXYDimensions(real, imag);
		final int length = real.length;
		Complex[] c = new Complex[length];
		for(int i = 0; i < length; ++i) {
			c[i] = new Complex(real[i], imag[i]);
		}
		return c;
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

	public static Complex mean(Complex[] a) {
		double realMean = 0.0;
		double imagMean = 0.0;

		for(int i = 0; i < a.length; ++i) {
			realMean += a[i].real();
			imagMean += a[i].imag();
		}
		return new Complex(realMean / a.length, imagMean / a.length);
	}

	public static Complex[] zeros(int length) {
		if(length < 0) {
			throw new IllegalArgumentException("dim must be greater than zero.");
		}
		Complex[] result = new Complex[length];
		for(int i = 0; i < length; ++i) {
			result[i] = new Complex();
		}
		return result;
	}

	public static Complex[] flatten(Complex[][] a) {
		int rows = a.length;
		int cols = a[0].length;
		Complex[] result = new Complex[rows * cols];
		for (int i = 0; i < rows; i++) {
			if(a[i].length != cols) {
				throw new IllegalArgumentException("All rows must have the same length.");
			}
			System.arraycopy(a[i], 0, result, i * cols, cols);
		}
		return result;
	}

	public static Complex[] add(Complex[] a, Complex b, int start, int end) {
		Complex[] result = new Complex[end - start];
		for(int i = start, k = 0; i < end; ++i, ++k) {
			result[k] = a[i].add(b);
		}
		return result;
	}

	public static Complex[] add(Complex[] a, Complex b) {
		return add(a, b, 0, a.length);
	}

	public static Complex[] reverse(Complex[] a) {
		final int length = a.length;
		Complex[] result = new Complex[length];
		for(int i = 0; i < length; ++i) {
			result[length - i - 1] = a[i];
		}
		return result;
	}

	public static Complex sum(Complex[] a) {
		final int length = a.length;
		Complex sum = new Complex();
		for(int i = 0; i < length; ++i) {
			sum.addEquals(a[i]);
		}
		return sum;
	}

	public static void multiplyInPlace(Complex[] a, double d) {
		for(int i = 0; i < a.length; ++i) {
			a[i].multiplyEquals(d);
		}
	}

	public static void divideInPlace(double d, Complex[] a) {
		for(int i = 0; i < a.length; ++i) {
			a[i].invertEquals();
			a[i].multiplyEquals(d);
		}
	}

	public static Complex[] divide(double d, Complex[] a) {
		Complex[] result = new Complex[a.length];
		for(int i = 0; i < a.length; ++i) {
			result[i] = a[i].invert().multiply(d);
		}
		return result;
	}

	public static Complex[] concatenate(Complex[] a, Complex[] b) {
		int aLength = a.length;
		int bLength = b.length;
		Complex[] result = new Complex[aLength + bLength];
		System.arraycopy(a, 0, result, 0, aLength);
		System.arraycopy(b, 0, result, aLength, bLength);
		return result;
	}

	public static Complex product(Complex[] a) {
		Complex prod = Complex.fromReal(1.0);
		for(int i = 0; i < a.length; ++i) {
			prod.multiplyEquals(a[i]);
		}
		return prod;
	}

	public static Complex[] multiply(Complex[] a, double d) {
		Complex[] result = new Complex[a.length];
		for(int i = 0; i < a.length; ++i) {
			result[i] = a[i].multiply(d);
		}
		return result;
	}

	public static Complex[] remove(Complex[] a, int index) {
		if(index < 0 || index >= a.length) {
			throw new IllegalArgumentException("index must be within the length of array a.");
		}
		Complex[] result = new Complex[a.length - 1];
		for(int i = 0, j = 0; i < result.length;) {
			if(j == index) {
				++j;
				continue;
			} else {
				result[i++] = a[j++];
			}
		}
		return result;
	}

	public static Complex[] remove(Complex[] a, Complex c) {
		int index = 0;
		for(int i = 0; i < a.length; ++i) {
			if(c.equals(a[i])) {
				index = i;
				break;
			}
		}
		return remove(a, index);
	}

	public static int[] argSort(Complex[] a) {
		Integer[] indexes = IntStream.range(0, a.length).boxed().toArray(Integer[]::new);
		Arrays.sort(indexes, Comparator.comparing(i -> a[i]));
		return Arrays.stream(indexes).mapToInt(Integer::intValue).toArray();
	}

	public static void main(String[] args) {
		Complex[] a = new Complex[]{
				Complex.fromReal(4.0),
				Complex.fromReal(3.0),
				Complex.fromReal(1.0),
				new Complex(2, -2),
				new Complex(2, 2),
				new Complex(2, -1),
				new Complex(2, 1),
				new Complex(2, -1),
				new Complex(2, 1),
				new Complex(1, 1),
				new Complex(1, -1)
		};

		Complex[] b = remove(a, 3);
		System.out.println(Arrays.toString(a));
		System.out.println(Arrays.toString(b));
		Complex[] c = remove(a, new Complex(2,2));
		System.out.println(Arrays.toString(c));

		System.out.println(Arrays.toString(argSort(a)));
	}

	public static Complex[] subtract(Complex[] a, Complex c) {
		Complex[] result = new Complex[a.length];
		for(int i = 0; i < a.length; ++i) {
			result[i] = a[i].subtract(c);
		}
		return result;
	}

	public static Complex[] subtract(Complex[] a, double d) {
		Complex[] result = new Complex[a.length];
		for(int i = 0; i < a.length; ++i) {
			result[i] = a[i].subtract(d);
		}
		return result;
	}

	public static Complex[] subtract(Complex c, Complex[] a) {
		Complex[] result = new Complex[a.length];
		for(int i = 0; i < a.length; ++i) {
			result[i] = c.subtract(a[i]);
		}
		return result;
	}
}
