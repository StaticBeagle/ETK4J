package com.wildbitsfoundry.etk4j.util;

import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;

import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.IntStream;

import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkXYDimensions;

public final class ComplexArrays {
    private ComplexArrays() {

    }

    /**
     * Deep copy of array.
     *
     * @param a The array to copy.
     * @return The deep copy i.e. new Complex values created from the real and imaginary part of the corresponding
     * Complex value in the input array.
     */
    public static Complex[] deepCopy(Complex[] a) {
        final int length = a.length;
        Complex[] result = new Complex[length];
        for (int i = 0; i < length; ++i) {
            result[i] = new Complex(a[i]);
        }
        return result;
    }

    /**
     * Convolve two arrays.
     *
     * @param a The left-hand array.
     * @param b The right-hand array.
     * @return The convolution of {@code a} and {@code b}.
     */
    public static Complex[] convolve(Complex[] a, Complex[] b) {
        Complex[] result = new Complex[a.length + b.length - 1];
        for (int i = 0; i < result.length; ++i) {
            result[i] = new Complex();
            for (int j = Math.max(0, i + 1 - b.length); j < Math.min(a.length, i + 1); ++j) {
                result[i].addEquals(a[j].multiply(b[i - j]));
            }
        }
        return result;
    }

    static int MOD = 998244353;

    /**
     * Real part of a Complex array element wise.
     *
     * @param a The input array.
     * @return An array of doubles containing the real part of each Complex value of the input array.
     */
    public static double[] real(Complex[] a) {
        double[] result = new double[a.length];
        for (int i = 0; i < a.length; ++i) {
            result[i] = a[i].real();
        }
        return result;
    }

    /**
     * Create a new Complex array from an array of real values.
     *
     * @param a The input array.
     * @return An array of Complex values with real part equal to the values of the input array and imaginary part
     * equal to zero.
     */
    public static Complex[] fromReal(double[] a) {
        Complex[] result = new Complex[a.length];
        for (int i = 0; i < a.length; ++i) {
            result[i] = Complex.fromReal(a[i]);
        }
        return result;
    }

    /**
     * Imaginary part of a Complex array element wise.
     *
     * @param a The input array.
     * @return An array of doubles containing the imaginary part of each Complex value of the input array.
     */
    public static double[] imag(Complex[] a) {
        double[] result = new double[a.length];
        for (int i = 0; i < a.length; ++i) {
            result[i] = a[i].imag();
        }
        return result;
    }

    /**
     * Create a new Complex array from an array of real values.
     *
     * @param a The input array.
     * @return An array of Complex values with imaginary part equal to the values of the input array and real part
     * equal to zero.
     */
    public static Complex[] fromImaginary(double[] a) {
        Complex[] result = new Complex[a.length];
        for (int i = 0; i < a.length; ++i) {
            result[i] = Complex.fromImaginary(a[i]);
        }
        return result;
    }

    /**
     * Convert and array representing the real part of an array of Complex number and another array representing the
     * imaginary part of an array of Complex into an array of Complex numbers.
     *
     * @param real The real part of the array of Complex to be created.
     * @param imag The imaginary part of the array of Complex to be created.
     * @return {@code {new Complex(real[i], imag[i], ..., new Complex(real[n], imag[n]}}.
     */
    public static Complex[] zip(double[] real, double[] imag) {
        checkXYDimensions(real, imag);
        final int length = real.length;
        Complex[] c = new Complex[length];
        for (int i = 0; i < length; ++i) {
            c[i] = new Complex(real[i], imag[i]);
        }
        return c;
    }

    /**
     * Mean of the array.
     *
     * @param a The input array.
     * @return The mean of the input array.
     */
    public static Complex mean(Complex[] a) {
        double realMean = 0.0;
        double imagMean = 0.0;

        for (Complex complex : a) {
            realMean += complex.real();
            imagMean += complex.imag();
        }
        return new Complex(realMean / a.length, imagMean / a.length);
    }

    /**
     * Creates and array of {@code n} 0 + j0.
     *
     * @param length The length of the array of ones.
     * @return An array containing {@code n} Complex(0, 0).
     */
    public static Complex[] zeros(int length) {
        if (length < 0) {
            throw new IllegalArgumentException("dim must be greater than zero.");
        }
        Complex[] result = new Complex[length];
        for (int i = 0; i < length; ++i) {
            result[i] = new Complex();
        }
        return result;
    }

    /**
     * Converts (flattens) a 2D array into a 1d array by copying
     * every row of the 2D array into the 1d array.
     * <pre>
     * let a = {{1, 2, 3}, {4, 5, 6}};
     * thus flatten(a) returns:
     * {1, 2, 3, 4, 5, 6}
     * </pre>
     *
     * @param a array to flatten
     * @return a new row-packed 1d array
     * @throws IllegalArgumentException if the input array a is jagged (i.e. not all rows have the same length)
     */
    public static Complex[] flatten(Complex[][] a) {
        int rows = a.length;
        int cols = a[0].length;
        Complex[] result = new Complex[rows * cols];
        for (int i = 0; i < rows; i++) {
            if (a[i].length != cols) {
                throw new IllegalArgumentException("All rows must have the same length.");
            }
            System.arraycopy(a[i], 0, result, i * cols, cols);
        }
        return result;
    }

    // region arithmetic operations

    /**
     * Add a Complex value to an array of Complex values, element wise.
     *
     * @param a The input array.
     * @param b The Complex value to be added.
     * @return {@code a + b}.
     **/
    public static Complex[] addElementWise(Complex[] a, Complex b) {
        Complex[] result = new Complex[a.length];
        for (int i = 0; i < a.length; ++i) {
            result[i] = a[i].add(b);
        }
        return result;
    }

    /**
     * Add two Complex arrays.
     * @param a The left-hard array.
     * @param b The right-hand array.
     * @return {@code a + b}.
     */
    public static Complex[] addElementWise(Complex[] a, Complex[] b) {
        int aLength = a.length;
        int bLength = b.length;
        if (aLength != bLength) {
            throw new IllegalArgumentException("Array dimensions must match.");
        }
        Complex[] result = new Complex[aLength];
        for (int i = 0; i < aLength; ++i) {
            result[i] = a[i].add(b[i]);
        }
        return result;
    }

    /**
     * Sum all the elements of the array element wise.
     * @param a The input array.
     * @return {@code sum(a)}.
     */
    public static Complex sum(Complex[] a) {
        final int length = a.length;
        Complex sum = new Complex();
        for (Complex complex : a) {
            sum.addEquals(complex);
        }
        return sum;
    }

    /**
     * Multiply and array of Complex values times a Complex value.
     * @param a The left-hand side array.
     * @param b The right-hand side array.
     * @return {@code a * b}
     */
    public static Complex[] multiplyElementWise(Complex[] a, Complex[] b) {
        Complex[] result = new Complex[a.length];
        for (int i = 0; i < a.length; ++i) {
            result[i] = a[i].multiply(b[i]);
        }
        return result;
    }

    /**
     * Multiply and array of Complex values times a Complex value.
     * @param a The left-hand side array.
     * @param b The Complex value to multiply with.
     * @return {@code a * b}
     */
    public static Complex[] multiplyElementWise(Complex[] a, Complex b) {
        Complex[] result = new Complex[a.length];
        for (int i = 0; i < a.length; ++i) {
            result[i] = a[i].multiply(b);
        }
        return result;
    }

    /**
     * Multiply and array of Complex values times a Complex value.
     * @param a The input array.
     * @param c The Complex value to multiply.
     * @return {@code a * b}
     */
    public static Complex[] multiplyElementWise(double[] a, Complex c) {
        Complex[] result = new Complex[a.length];
        for (int i = 0; i < a.length; ++i) {
            result[i] = c.multiply(a[i]);
        }
        return result;
    }

    /**
     * Multiply and array of Complex values times a real value.
     * @param a The input array.
     * @param d The real value to multiply.
     * @return {@code a * b}
     */
    public static Complex[] multiplyElementWise(Complex[] a, double d) {
        Complex[] result = ComplexArrays.deepCopy(a);
        multiplyElementWiseInPlace(result, d);
        return result;
    }

    /**
     * Multiply and array of Complex values times a real value in place. The result of the multiplication is stored in
     * the input array a.
     * @param a The input array.
     * @param d The real value to multiply.
     */
    public static void multiplyElementWiseInPlace(Complex[] a, double d) {
        for (Complex complex : a) {
            complex.multiplyEquals(d);
        }
    }

    /**
     * Divides a real number by an array of Complex number element wise.
     * @param d The left-hand double value.
     * @param a The right-hand array of Complex values.
     * @return {@code d / a}.
     */
    public static Complex[] divideElementWise(double d, Complex[] a) {
        Complex[] result = new Complex[a.length];
        for (int i = 0; i < a.length; ++i) {
            result[i] = a[i].invert().multiply(d);
        }
        return result;
    }

    /**
     * Divides a real number by an array of Complex number element wise.
     * @param a The left-hand array of Complex values.
     * @param d The right-hand double value.
     * @return {@code a / d}.
     */
    public static Complex[] divideElementWise(Complex[] a, double d) {
        Complex[] result = new Complex[a.length];
        double dp = 1.0 / d;
        for (int i = 0; i < a.length; ++i) {
            result[i] = a[i].multiply(dp);
        }
        return result;
    }

    /**
     * Array concatenation.
     * @param a The left-hand array.
     * @param b The right-hand array.
     * @return {@code {a, b}}
     */
    public static Complex[] concatenate(Complex[] a, Complex[] b) {
        int aLength = a.length;
        int bLength = b.length;
        Complex[] result = new Complex[aLength + bLength];
        System.arraycopy(a, 0, result, 0, aLength);
        System.arraycopy(b, 0, result, aLength, bLength);
        return result;
    }

    /**
     * Array product.
     * @param a The input array.
     * @return The product of multiplying all the elements in the array.
     */
    public static Complex product(Complex[] a) {
        Complex prod = Complex.fromReal(1.0);
        for (Complex complex : a) {
            prod.multiplyEquals(complex);
        }
        return prod;
    }

    /**
     * Subtract a Complex value from an array of Complex values element wise.
     *
     * @param a The left-hand array.
     * @param b The right-hand array to be subtracted.
     * @return {@code a - b}.
     **/
    public static Complex[] subtractElementWise(Complex[] a, Complex[] b) {
        checkXYDimensions(a, b);
        Complex[] result = new Complex[a.length];
        for (int i = 0; i < a.length; ++i) {
            result[i] = a[i].subtract(b[i]);
        }
        return result;
    }

    /**
     * Subtract a Complex value from an array of Complex values element wise.
     *
     * @param a The left-hand array.
     * @param c The right-hand Complex value to be subtracted.
     * @return {@code a - b}.
     **/
    public static Complex[] subtractElementWise(Complex[] a, Complex c) {
        Complex[] result = new Complex[a.length];
        for (int i = 0; i < a.length; ++i) {
            result[i] = a[i].subtract(c);
        }
        return result;
    }

    /**
     * Subtract a real value from an array of Complex values element wise.
     *
     * @param a The left-hand array.
     * @param d The right-hand real value to be subtracted.
     * @return {@code a - b}.
     **/
    public static Complex[] subtractElementWise(Complex[] a, double d) {
        Complex[] result = new Complex[a.length];
        for (int i = 0; i < a.length; ++i) {
            result[i] = a[i].subtract(d);
        }
        return result;
    }

    /**
     * Subtract an array of Complex values from a Complex number element wise.
     *
     * @param c The Complex number.
     * @param a The right-hand array to be subtracted.
     * @return {@code c - a}.
     **/
    public static Complex[] subtractElementWise(Complex c, Complex[] a) {
        Complex[] result = new Complex[a.length];
        for (int i = 0; i < a.length; ++i) {
            result[i] = c.subtract(a[i]);
        }
        return result;
    }

    /**
     * Norm 2 of the Complex array.
     *
     * @param x The Complex array.
     * @return (&Sigma;|x<sub>i</sub>|)<sup>½</sup>
     **/
    public static double norm(Complex[] x) {
        double nrm = 0.0;
        // Compute 2-norm without under/overflow.
        for (Complex complex : x) {
            nrm = MathETK.hypot(nrm, complex.abs());
        }
        return nrm;
    }

    public static double normFro(Complex[] a) {
        int i;
        double fac, nrm, scale;

        int n = a.length;

        scale = 0.0;
        for (i = 0; i < n; i++) {
            scale = Math.max(scale,
                    Math.abs(a[i].real()) + Math.abs(a[i].imag()));
        }
        if (scale == 0) {
            return 0.0;
        }
        if (scale < 1) {
            scale = scale * 1.0e20;
        }
        scale = 1 / scale;
        nrm = 0;

        for (i = 0; i < n; i++) {
            fac = scale * a[i].real();
            nrm = nrm + fac * fac;
            fac = scale * a[i].imag();
            nrm = nrm + fac * fac;
        }

        return Math.sqrt(nrm) / scale;
    }
    
    public static Complex[][] zeros(int rows, int cols) {
        Complex[][] zeros = new Complex[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                zeros[i][j] = new Complex();
            }
        }
        return zeros;
    }
}
