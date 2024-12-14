package com.wildbitsfoundry.etk4j.util;

import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.IntStream;

import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.interpolation.CubicSpline;

/**
 * The {@code DoubleArrays} utility class provides methods to manipulate arrays of native {@code double} values.
 */
public final class DoubleArrays {

    private DoubleArrays() {
    }

    /***
     * Creates n linearly spaced samples between x0 and x1.
     *
     * @param x0
     *            starting point
     * @param x1
     *            end point
     * @param n
     *            number of samples
     * @return Array of n equally spaced samples between [x0, x1]
     */
    public static double[] linSpace(double x0, double x1, int n) {
        double[] result = new double[n];
        int i = 0;
        double delta = (x1 - x0) / --n;
        while (i < n) {
            result[i] = x0 + i * delta;
            ++i;
        }
        result[i] = x1;
        return result;
    }

    /***
     * Creates n linearly spaced samples between x0 and x1
     *
     * @param x0
     *            starting point
     * @param x1
     *            end point
     * @param step
     *            step size
     * @return Array of n equally spaced samples between [x0, x1]
     */
    public static double[] linSteps(double x0, double x1, double step) {
        final int n = (int) Math.ceil((x1 - x0) / step);
        double[] result = new double[n + 1];
        int i = 0;
        while (i < n) {
            result[i] = x0 + i * step;
            ++i;
        }
        result[i] = x1;
        return result;
    }

    /***
     * Creates n linearly spaced samples between x0 and x1
     *
     * @param x0
     *            starting point
     * @param x1
     *            end point
     * @return Array of n equally spaced samples between [x0, x1]. Step size is
     *         assumed to be 1
     */
    public static double[] linSteps(double x0, double x1) {
        return linSteps(x0, x1, 1.0);
    }

    /***
     * Creates n logarithmically spaced samples between decades x0 and x1
     *
     * @param x0
     *            starting decade
     * @param x1
     *            end decade
     * @param n
     *            number of samples
     * @return Array of n logarithmically spaced samples between [x0, x1]
     */
    public static double[] logSpace(int x0, int x1, int n) {
        double[] result = new double[n];
        int i = 0;
        double delta = (double) (x1 - x0) / --n;
        while (i < n) {
            result[i] = Math.pow(10.0, x0 + i * delta);
            ++i;
        }
        result[i] = Math.pow(10.0, x1);
        return result;
    }

    /**
     * Subtract two arrays element wise.
     * @param a The left-hand array.
     * @param b The right-hand array.
     * @return {@code a - b}.
     */
    public static double[] subtractElementWise(double[] a, double[] b) {
        double[] result = Arrays.copyOf(a, a.length);
        subtractElementWiseInPlace(result, b);
        return result;
    }

    /**
     * Subtract an array element wise from a scalar.
     * @param a The scalar argument.
     * @param b The array to subtract.
     * @return {@code a - b}.
     */
    public static double[] subtractElementWise(double a, double[] b) {
        final int length = b.length;
        double[] result = new double[length];
        for (int i = 0; i < length; ++i) {
            result[i] = a - b[i];
        }
        return result;
    }

    /**
     * Subtract a scalar element wise from array a.
     * @param a The left-hand array.
     * @param b The scalar to subtract.
     * @return {@code a - b}.
     */
    public static double[] subtractElementWise(double[] a, double b) {
        final int length = a.length;
        double[] result = new double[length];
        for (int i = 0; i < length; ++i) {
            result[i] = a[i] - b;
        }
        return result;
    }

    /**
     * Subtract two arrays element wise in place. The result of the subtraction is stored in array {@code a}.
     * @param a The left-hand array.
     * @param b The right-hand array.
     */
    public static void subtractElementWiseInPlace(double[] a, double[] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b dimensions must match.");
        }
        for (int i = 0; i < a.length; ++i) {
            a[i] = a[i] - b[i];
        }
    }


    /**
     * Convolve two arrays.
     * @param a The left-hand array.
     * @param b The right-hand array.
     * @return The convolution of {@code a} and {@code b}.
     */
    public static double[] convolve(double[] a, double[] b) {
        final int n = a.length;
        final int m = b.length;

        double[] result = new double[n + m - 1];
        for(int i = 0; i < n; ++i)
        {
            for(int j = 0; j < m; ++j)
            {
                result[i + j] += a[i] * b[j];
            }
        }
        return result;
    }

    /**
     * Add an array and a scalar element wise.
     * @param a The array to multiply.
     * @param b The scalar.
     * @return {@code a + b}.
     */
    public static double[] addElementWise(double[] a, double b) {
        double[] result = Arrays.copyOf(a, a.length);
        addElementWiseInPlace(result, b);
        return result;
    }

    /**
     * Add two arrays element wise.
     * @param a The left-hand array.
     * @param b The right-hand array.
     * @return {@code a + b}.
     */
    public static double[] addElementWise(double[] a, double[] b) {
        double[] result = Arrays.copyOf(a, a.length);
        addElementWiseInPlace(result, b);
        return result;
    }

    /**
     * Add an array and a scalar element wise in place. The result of the addition is stored in array a.
     * @param a The array to multiply.
     * @param b The scalar.
     */
    public static void addElementWiseInPlace(double[] a, double b) {
        for (int i = 0; i < a.length; ++i) {
            a[i] += b;
        }
    }

    /**
     * Add two arrays element wise in place. The result of the addition is stored in array a.
     * @param a The left-hand array.
     * @param b The right-hand array.
     */
    public static void addElementWiseInPlace(double[] a, double[] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b dimensions must match.");
        }
        for (int i = 0; i < a.length; ++i) {
            a[i] += b[i];
        }
    }


    /**
     * Multiply two and array a scalar element wise.
     * @param a The array to multiply.
     * @param d The scalar to multiply.
     * @return {@code a * b}.
     */
    public static double[] multiplyElementWise(double[] a, double d) {
        double[] result = Arrays.copyOf(a, a.length);
        multiplyElementWiseInPlace(result, d);
        return result;
    }

    /**
     * Multiply two and array a scalar element wise in place. The result of the multiplication is stored in array a.
     * @param a The array to multiply.
     * @param d The scalar to multiply.
     */
    public static void multiplyElementWiseInPlace(double[] a, double d) {
        final int length = a.length;
        for (int i = 0; i < length; ++i) {
            a[i] *= d;
        }
    }

    /**
     * Multiply two arrays element wise.
     * @param a The left-hand array.
     * @param b The right-hand array.
     * @return {@code a * b}.
     */
    public static double[] multiplyElementWise(double[] a, double[] b) {
        double[] result = Arrays.copyOf(a, a.length);
        multiplyElementWiseInPlace(result, b);
        return result;
    }

    /**
     * Multiply two arrays element wise in place. The result of the multiplication is stored in array a.
     * @param a The left-hand array.
     * @param b The right-hand array.
     */
    public static void multiplyElementWiseInPlace(double[] a, double[] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b dimensions must match.");
        }
        final int length = a.length;
        for (int i = 0; i < length; ++i) {
            a[i] *= b[i];
        }
    }
//TODO test
    /**
     * Array division element wise.
     * @param a The left-hand array.
     * @param d The divisor.
     * @return {@code a / b}.
     */
    public static double[] divideElementWise(double[] a, double d) {
        double[] result = Arrays.copyOf(a, a.length);
        divideElementWiseInPlace(result, d);
        return result;
    }

    /**
     * Array division element wise.
     * @param a The left-hand array.
     * @param b The right-hand array.
     * @return {@code a / b}.
     */
    public static double[] divideElementWise(double[] a, double[] b) {
        double[] result = Arrays.copyOf(a, a.length);
        divideElementWiseInPlace(result, b);
        return result;
    }
// TODO test
    /**
     * Array division element wise in place. The result of the division is stored in array a.
     * @param a The left-hand array.
     * @param d The divisor.
     */
    public static void divideElementWiseInPlace(double[] a, double d) {
        multiplyElementWiseInPlace(a, 1 / d);
    }

    /**
     * Array division element wise in place. The result of the division is stored in array a.
     * @param a The left-hand array.
     * @param b The right-hand array.
     */
    public static void divideElementWiseInPlace(double[] a, double[] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b dimensions must match.");
        }

        final int length = a.length;
        for (int i = 0; i < length; ++i) {
            a[i] /= b[i];
        }
    }

    /**
     * Maximum value in the array.
     * @param a The array whose maximum has to be found.
     * @return {@code max(a)}
     */
    public static double max(double... a) {
        double max = a[0];
        for (int i = 1; i < a.length; ++i) {
            if (a[i] > max) {
                max = a[i];
            }
        }
        return max;
    }

    /**
     * Minimum value in the array.
     * @param a The array whose minimum has to be found.
     * @return {@code min(a)}
     */
    public static double min(double... a) {
        double min = a[0];
        for (int i = 1; i < a.length; ++i) {
            if (a[i] < min) {
                min = a[i];
            }
        }
        return min;
    }

    /**
     * Index of minimum value.
     * @param a The array whose index of minimum value has to be found.
     * @return The index at which the minimum of array {@code a} occurs.
     */
    public static int argMin(double[] a) {
        int index = 0;
        double min = a[0];
        for (int i = 1; i < a.length; ++i) {
            if (a[i] < min) {
                min = a[i];
                index = i;
            }
        }
        return index;
    }

    /**
     * Returns the indices that would sort an array.
     * @param a The input array.
     * @return The indexes that would sort the array {@code a}.
     */
    public static int[] argSort(double[] a) {
        Integer[] indexes = IntStream.range(0, a.length).boxed().toArray(Integer[]::new);
        Arrays.sort(indexes, Comparator.comparing(i -> a[i]));
        return Arrays.stream(indexes).mapToInt(Integer::intValue).toArray();
    }

    /**
     * Norm one of array.
     * @param a The array whose norm needs to be found.
     * @return {@code sum(abs(a))}.
     */
    public static double norm1(double[] a) {
        double norm = 0.0;
        for (int i = 0; i < a.length; ++i) {
            norm += Math.abs(a[i]);
        }
        return norm;
    }

    /**
     * Norm two of array.
     * @param a The array whose norm needs to be found.
     * @return {@code sum(abs(a)<sup>2</sup>)<sup>1/2</sup>}.
     */
    public static double norm2(double[] a) {
        double max = normInf(a);
        if (max == 0.0) {
            return 0.0;
        }
        double tmp = Math.pow(1.0 / max, 2);
        double norm = 0.0;
        for (int i = 0; i < a.length; ++i) {
            norm += a[i] * a[i] * tmp;
        }
        return Math.sqrt(norm) * max;
    }

    /**
     * Infinite norm of array.
     * @param a The array whose norm needs to be found.
     * @return {@code max(abs(a))}.
     */
    public static double normInf(double[] a) {
        double max = Math.abs(a[0]);
        for (int i = 1; i < a.length; ++i) {
            double abs = Math.abs(a[i]);
            if (abs > max) {
                max = abs;
            }
        }
        return max;
    }

    /**
     * Negative Infinite norm of array.
     * @param a The array whose norm needs to be found.
     * @return {@code min(abs(a))}.
     */
    public static double normNegInf(double[] a) {
        double min = Math.abs(a[0]);
        for (int i = 1; i < a.length; ++i) {
            double abs = Math.abs(a[i]);
            if (abs < min) {
                min = abs;
            }
        }
        return min;
    }

    /***
     * Euclidean distance between two arrays.
     *
     * @param a The left-hand array.
     * @param b The right-hand array.
     * @return {@code sqrt(sum((a[i] - b[i])<sup>2</sup>))}
     */
    public static double distance(double[] a, double[] b) {
        double norm = 0.0;
        for (int i = 0; i < a.length; ++i) {
            norm = MathETK.hypot(norm, a[i] - b[i]);
        }
        return norm;
    }

    /***
     * Kronecker array product.
     * @param a input array
     * @param b input array
     * @return an array formed by taking all possible products between the elements of a and b
     */
    public static double[] kronecker(final double[] a, final double[] b) {
        int aLength = a.length;
        int bLength = b.length;
        double[] c = new double[aLength * bLength];
        if (a.length == 1) {
            return DoubleArrays.multiplyElementWise(b, a[0]);
        }
        if (b.length == 1) {
            return DoubleArrays.multiplyElementWise(a, b[0]);
        }
        for (int i = 0; i < aLength; i++) {
            for (int j = 0; j < bLength; j++) {
                c[j + i * (bLength)] = a[i] * b[j];
            }
        }
        return c;
    }

    /**
     * Array concatenation.
     * @param a The left-hand array.
     * @param b The right-hand array.
     * @return {@code {a, b}}
     */
    public static double[] concatenate(double[] a, double[] b) {
        int aLength = a.length;
        int bLength = b.length;
        double[] result = new double[aLength + bLength];
        System.arraycopy(a, 0, result, 0, aLength);
        System.arraycopy(b, 0, result, aLength, bLength);
        return result;
    }

    /**
     * Concatenate series of arrays.
     * @param a The left-hand array.
     * @param rest The rest of the arrays to concatenate.
     * @return {@code {a, rest[0], rest[1], ... , rest[n}}
     */
    public static double[] concatenateAll(double[] a, double[]... rest) {
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

    /**
     * Repeat an array.
     * @param a The array to repeat.
     * @param n The number of times to repeat the array.
     * @return The input array repeated {@code n} times.
     */
    public static double[] repeat(double[] a, int n) {
        final int length = a.length;
        double[] result = new double[length * n];
        for (int i = 0; i < n; ++i) {
            System.arraycopy(a, 0, result, i * length, length);
        }
        return result;
    }

    /**
     * Reverse and array.
     * @param a The array to reverse.
     * @return The input array reversed.
     */
    public static double[] reverse(double[] a) {
        final int length = a.length;
        double[] result = new double[a.length];
        for (int i = 0, j = length - 1; i <= length / 2; ++i, --j) {
            result[i] = a[j];
            result[j] = a[i];
        }
        return result;
    }

    /**
     * The sum of all the elements in the array.
     * @param a The array whose sum needs to be found.
     * @return {@code sum(a)}.
     */
    public static double sum(double[] a) {
        final int length = a.length;
        double result = 0;
        for (int i = 0; i < length; ++i) {
            result += a[i];
        }
        return result;
    }

    /**
     * The sum of the square of all the elements in the array.
     * @param a The array whose sum of all squares needs to be found.
     * @return {@code sum(a)}
     */
    public static double sumSquares(double[] a) {
        final int length = a.length;
        double result = 0;
        for (int i = 0; i < length; ++i) {
            result += a[i] * a[i];
        }
        return result;
    }

    /**
     * Cumulative sum.
     * @param a The array whose cumulative sum needs to found.
     * @return The cumulative sum of a.
     */
    public static double[] cumulativeSum(double[] a) {
        final int length = a.length;
        double[] result = new double[length];
        result[0] = a[0];
        for (int i = 1; i < length; ++i) {
            result[i] = result[i - 1] + a[i];
        }
        return result;
    }

    /**
     * Root Mean Square.
     * @param a The array whose RMS value needs to be found.
     * @return {@code rms(a)}.
     */
    public static double rms(double[] a) {
        final int N = a.length;
        return norm2(a) / Math.sqrt(N);
    }

    /**
     * Is array ascending.
     * @param a The input array.
     * @return {@code true} if all the elements in the array are sorted in ascending order, {@code false} otherwise.
     */
    public static boolean isAscending(double[] a) {
        final int n = a.length;
        for (int i = 0; i < n - 1; ++i) {
            if (a[i] > a[i + 1]) {
                return false;
            }
        }
        return true;
    }


    /**
     * Converts (flattens) a 2D array into a 1d array by copying
     * every row of the 2D array into the 1d array.
     * <pre>
     * let a = {{1, 2, 3}, {4, 5, 6}};
     * thus flatten(a) returns:
     * {1, 2, 3, 4, 5, 6}
     *</pre>
     * @param a array to flatten
     * @return a new row-packed 1d array
     * @throws IllegalArgumentException if the input array a is jagged (i.e. not all rows have the same length)
     */
    public static double[] flatten(double[][] a) {
        int rows = a.length;
        int cols = a[0].length;
        double[] result = new double[rows * cols];
        for (int i = 0; i < rows; i++) {
            if(a[i].length != cols) {
                throw new IllegalArgumentException("All rows must have the same length.");
            }
            System.arraycopy(a[i], 0, result, i * cols, cols);
        }
        return result;
    }

    /**
     * Mean of the array.
     * @param y The input array.
     * @return The mean of the input array.
     */
    public static double mean(double[] y) {
        final int n = y.length;
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            sum += y[i];
        }
        return sum / n;
    }

    /**
     * Horner's method of evaluating a polynomial.
     * @param coefficients The coefficients of the polynomial.
     * @param x Argument at which to evaluate the polynomial.
     * @return The polynomial evaluated at {@code x}.
     */
    public static double horner(double[] coefficients, double x) {
        double result = 0.0;
        for (int j = 0; j < coefficients.length; ++j) {
            result = result * x + coefficients[j];
        }
        return result;
    }

    /**
     * Differences between adjacent elements in the array.
     * @param a The input array.
     * @return The difference between adjacent elements in the array.
     */
    public static double[] difference(double[] a) {
        double[] result = new double[a.length - 1];
        for(int i = 0; i < result.length; ++i) {
            result[i] = a[i + 1] - a[i];
        }
        return result;
    }

    /**
     * Numerical gradient of an array. The spacing between the elements is assumed to be one.
     * @param a The input array.
     * @return The numerical gradient of the array.
     */
    public static double[] gradient(double[] a) {
        double[] grad = new double[a.length];
        final int end = a.length - 1;
        grad[0] = a[1] - a[0]; // one-sided difference
        grad[end] = a[end] - a[end - 1]; // one-sided difference

        // central difference
        for (int i = 1; i < end; ++i) {
            grad[i] = 0.5 * (a[i + 1] - a[i - 1]);
        }
        return grad;
    }

    /**
     * Numerical gradient of an array.
     * @param a The input array.
     * @param h The spacing between each value in {@code a}.
     * @return The numerical gradient of the array.
     */
    public static double[] gradient(double[] a, double[] h) {
        if(a.length != h.length) {
            throw new IllegalArgumentException("a and h dimensions must match.");
        }

        double[] grad = new double[a.length];
        final int end = a.length - 1;
        grad[0] = (a[1] - a[0]) / (h[1] - h[0]); // one-sided difference
        grad[end] = (a[end] - a[end - 1]) / (h[end] - h[end - 1]); // one-sided difference

        // central difference
        for (int i = 1; i < end; ++i) {
            grad[i] = (a[i + 1] - a[i - 1]) / (h[i + 1] - h[i - 1]);
        }
        return grad;
    }

    /**
     * Outer product of two arrays.
     * @param a The first array.
     * @param b The second array.
     * @return {@code output[i, j] = a[i] * b[j]}.
     */
    public static double[][] outer(double[] a, double[] b) {
        int M = a.length;
        int N = b.length;
        double[][] result = new double[M][N];
        for(int i = 0; i < M; ++i) {
            for(int j = 0; j < N; ++j) {
                result[i][j] = a[i] * b[j];
            }
        }
        return result;
    }

    /**
     * Dot product between two arrays.
     * @param a The first array.
     * @param b The second array.
     * @return The dot (inner) product of arrays {@code a} and {@code b}.
     */
    public static double dot(double[] a, double[] b) {
        if(a.length != b.length) {
            throw new IllegalArgumentException("Both arrays must be of the same length.");
        }
        double result = 0.0;
        for(int i = 0; i < a.length; ++i) {
            result += a[i] * b[i];
        }
        return result;
    }

    /**
     * Dot product between multiple arrays.
     * @param a Array of arrays.
     * @param b The second array.
     * @return The dot (inner) product of arrays {@code a[i]} and {@code b}.
     */
    public static double[] dot(double[][] a, double[] b) {
        double[] result = new double[a.length];
        for(int i = 0; i < a.length; i++) {
            result[i] = dot(a[i], b);
        }
        return result;
    }

    /**
     * Array product.
     * @param a The input array.
     * @return The product of multiplying all the elements in the array.
     */
    public static double product(double[] a) {
        double prod = 1.0;
        for(int i = 0; i < a.length; ++i) {
            prod *= a[i];
        }
        return prod;
    }

    /**
     * Checks if all elements of the array are close numerically to a given target.
     *
     * @param a The input array.
     * @param target The target value.
     * @return {@code true} if all the values in the array are close to the given target, {@code false} otherwise.
     * @see MathETK#isClose(double, double) isClose(double, double).
     */
    public static boolean allClose(double[] a, double target) {
        for(int i = 0; i < a.length; ++i) {
            if(!MathETK.isClose(a[i], target)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Checks if all elements of the array are close numerically to a given target.
     *
     * @param a The input array.
     * @param target The target value.
     * @param absTol The absolute tolerance.
     * @param relTol The relative tolerance.
     * @return {@code true} if all the values in the array are close to the given target, {@code false} otherwise.
     * @see MathETK#isClose(double, double, double, double) isClose(double, double, double, double).
     */
    public static boolean allClose(double[] a, double target, double absTol, double relTol) {
        for(int i = 0; i < a.length; ++i) {
            if(!MathETK.isClose(a[i], target, absTol, relTol)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Checks if all elements of the array are close numerically to a given target.
     *
     * @param a The input array.
     * @param target The target value.
     * @param absTol The absolute tolerance.
     * @return {@code true} if all the values in the array are close to the given target, {@code false} otherwise.
     * @see MathETK#isClose(double, double, double) isClose(double, double, double).
     */
    public static boolean allClose(double[] a, double target, double absTol) {
        for(int i = 0; i < a.length; ++i) {
            if(!MathETK.isClose(a[i], target, absTol)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Creates and array of {@code n} ones.
     * @param n The length of the array of ones.
     * @return An array containing {@code n} double ones.
     */
    public static double[] ones(int n) {
        double[] result = new double[n];
        Arrays.fill(result, 1.0);
        return result;
    }

    /**
     * Transpose a 2D array.
     * @param a The array to transpose.
     * @return The transposed array.
     * @see <a href="https://stackoverflow.com/a/17634025/6383857">Transpose 2D array</a>
     */
    public static double[][] transpose(double[][] a) {
        // empty or unset array, nothing do to here
        if (a == null || a.length == 0)
            return a;

        int width = a.length;
        int height = a[0].length;

        double[][] result = new double[height][width];

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                result[y][x] = a[x][y];
            }
        }
        return result;
    }

    /***
     * Calculates the absolute value element wise for all the elements in the array
     * @param a The array
     * @return {@code abs(a)}
     */
    public static double[] abs(double[] a) {
        double[] result = new double[a.length];
        for(int i = 0; i < a.length; i++) {
            result[i] = Math.abs(a[i]);
        }
        return result;
    }

    public static double[] max(double[] a, double[] b) {
        if(a.length != b.length) {
            throw new IllegalArgumentException("The length of the two arrays must be the same");
        }
        double[] result = new double[a.length];
        for(int i = 0; i < a.length; i++) {
            result[i] = Math.max(a[i], b[i]);
        }
        return result;
    }
}
