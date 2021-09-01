package com.wildbitsfoundry.etk4j.util;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.interpolation.CubicSpline;

public final class NumArrays {

    private NumArrays() {
    }

    /***
     * Creates n linearly spaced samples between x0 and x1
     *
     * @param x0
     *            starting point
     * @param x1
     *            end point
     * @param n
     *            number of samples
     * @return Array of n equally spaced samples between [x0, x1]
     */
    public static double[] linspace(double x0, double x1, int n) {
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
     * @param step
     *            step size
     * @param x1
     *            end point
     * @return Array of n equally spaced samples between [x0, x1]
     */
    public static double[] linsteps(double x0, double step, double x1) {
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
    public static double[] linsteps(double x0, double x1) {
        return linsteps(x0, 1.0, x1);
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
    public static double[] logspace(int x0, int x1, int n) {
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

    public static double[] add(double[] a, double[] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b dimensions must match");
        }
        final int length = a.length;
        double[] result = new double[length];
        for (int i = 0; i < length; ++i) {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    public static double[] subtract(double[] a, double[] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b dimensions must match");
        }
        final int length = a.length;
        double[] result = new double[length];
        for (int i = 0; i < length; ++i) {
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

    public static double[] addElementWise(double[] a, double b) {
        double[] result = Arrays.copyOf(a, a.length);
        addElementWiseInPlace(result, b);
        return result;
    }

    public static void addElementWiseInPlace(double[] a, double b) {
        for (int i = 0; i < a.length; ++i) {
            a[i] += b;
        }
    }

    public static double[] multiplyElementWise(double[] a, double[] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b dimensions must match");
        }

        double[] result = Arrays.copyOf(a, a.length);
        multiplyElementWiseInPlace(result, b);
        return result;
    }

    public static void multiplyElementWiseInPlace(double[] a, double[] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b dimensions must match");
        }

        final int length = a.length;
        for (int i = 0; i < length; ++i) {
            a[i] *= b[i];
        }
    }

    public static double[] divideElementWise(double[] a, double[] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b dimensions must match");
        }

        double[] result = Arrays.copyOf(a, a.length);
        divideElementWiseInPlace(result, b);
        return result;
    }

    public static void divideElementWiseInPlace(double[] a, double[] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b dimensions must match");
        }

        final int length = a.length;
        for (int i = 0; i < length; ++i) {
            a[i] /= b[i];
        }
    }

    public static double normFast(double[] a) {
        double norm = 0.0;
        for (int i = 0; i < a.length; ++i) {
            norm += a[i] * a[i];
        }
        return Math.sqrt(norm);
    }

    public static double max(double... a) {
        double max = a[0];
        for (int i = 1; i < a.length; ++i) {
            if (a[i] > max) {
                max = a[i];
            }
        }
        return max;
    }

    public static double min(double... a) {
        double min = a[0];
        for (int i = 1; i < a.length; ++i) {
            if (a[i] < min) {
                min = a[i];
            }
        }
        return min;
    }

    public static double norm1(double[] a) {
        double norm = 0.0;
        for (int i = 0; i < a.length; ++i) {
            norm += Math.abs(a[i]);
        }
        return norm;
    }

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

    public static double normInf(Double[] a) {
        double max = Math.abs(a[0].doubleValue());
        for (int i = 1; i < a.length; ++i) {
            double abs = Math.abs(a[i].doubleValue());
            if (abs > max) {
                max = abs;
            }
        }
        return max;
    }

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
     * Euclidean distance between two arrays
     *
     * @param a
     *            - first array
     * @param b
     *            - second array
     * @return
     *
     *         <pre>
     *         sqrt(sum((a[i] - b[i])<sup>2</sup>))
     *         </pre>
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
        for (int i = 0; i < n; ++i) {
            System.arraycopy(a, 0, result, i * length, length);
        }
        return result;
    }

    public static double[] reverse(double[] a) {
        final int length = a.length;
        double[] result = new double[a.length];
        for (int i = length - 1, j = 0; i >= 0; --i, ++j) {
            result[j] = a[i];
        }
        return result;
    }

    public static double sum(double[] a) {
        final int length = a.length;
        double result = 0;
        for (int i = 0; i < length; ++i) {
            result += a[i];
        }
        return result;
    }

    public static double sumSquare(double[] a) {
        final int length = a.length;
        double result = 0;
        for (int i = 0; i < length; ++i) {
            result += a[i] * a[i];
        }
        return result;
    }

    public static double[] cummulativeSum(double[] a) {
        final int length = a.length;
        double[] result = new double[length];
        result[0] = a[0];
        for (int i = 1; i < length; ++i) {
            result[i] = result[i - 1] + a[i];
        }
        return result;
    }

    public static double rms(double[] a) {
        final int N = a.length;
        return norm2(a) / Math.sqrt(N);
    }

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
     * Converts (flattens) a 2d array into a 1d array by copying
     * every row of the 2d array into the 1d array.
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

    public static double mean(double[] y) {
        final int n = y.length;
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            sum += y[i];
        }
        return sum / n;
    }

    public static Double[] box(double[] a) {
        final int length = a.length;
        Double[] boxed = new Double[length];
        for (int i = 0; i < length; ++i) {
            boxed[i] = a[i];
        }
        return boxed;
    }

    public static double[] unbox(Double[] a) {
        final int length = a.length;
        double[] unboxed = new double[length];
        for (int i = 0; i < length; ++i) {
            unboxed[i] = a[i].doubleValue();
        }
        return unboxed;
    }

    public static double[] scaleRange(double[] a, double curMin, double curMax, double newMin, double newMax) {
        double[] result = new double[a.length];
        scaleRangeInPlace(result, curMax, curMin, newMax, newMin);
        return result;
    }

    public static void scaleRangeInPlace(double[] a, double curMin, double curMax, double newMin, double newMax) {
        final int length = a.length;
        double delta = (newMax - newMin) / (curMax - curMin);
        for (int i = 0; i < length; ++i) {
            a[i] = delta * (a[i] - curMin) + curMin;
        }
    }

    public static double horner(double[] coefs, double x) {
        double result = 0.0;
        for (int j = 0; j < coefs.length; ++j) {
            result = result * x + coefs[j];
        }
        return result;
    }

    public static boolean equals(double[] a, double d) {
        for (double val : a) {
            if (d != val) {
                return false;
            }
        }
        return true;
    }

    public static double[] gradient(double[] a) {
        // check bounds

        double[] grad = new double[a.length];
        final int end = a.length - 1;
        grad[0] = a[1] - a[0]; // once sided difference
        grad[end] = a[end] - a[end - 1]; // once sided difference

        // central difference
        for (int i = 1; i < end; ++i) {
            grad[i] = 0.5 * (a[i + 1] - a[i - 1]);
        }
        return grad;
    }

    public static double[] gradient(double[] a, double[] h) {
        // check bounds
        // check that arrays have the same size

        double[] grad = new double[a.length];
        final int end = a.length - 1;
        grad[0] = (a[1] - a[0]) / (h[1] - h[0]); // once sided difference
        grad[end] = (a[end] - a[end - 1]) / (h[end] - h[end - 1]); // once sided difference

        // central difference
        for (int i = 1; i < end; ++i) {
            grad[i] = (a[i + 1] - a[i - 1]) / (h[i + 1] - h[i - 1]);
        }
        return grad;
    }

    public static double[] expand(double[] a, int newSize) {
        return resize(a, newSize, (double) a.length / (newSize + 1));
    }

    public static double[] shrink(double[] a, int newSize) {
        return resize(a, newSize, (double) a.length / (newSize - 1));
    }

    private static double[] resize(double[] a, int newSize, double step) {
        final int length = a.length;

        double[] x = new double[length];
        for (int i = 0; i < length; ++i) {
            x[i] = i;
        }
        CubicSpline sp = CubicSpline.newNotAKnotSpline(x, a);

        double[] result = new double[newSize];
        double xi = x[0];
        for (int i = 0; i < newSize - 1; ++i) {
            result[i] = sp.evaluateAt(xi);
            xi += step;
        }
        result[newSize - 1] = a[length - 1];
        return result;
    }

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

    public static double dot(double[] a, double[] b) {
        if(a.length != b.length) {
            throw new IllegalArgumentException("Both arrays must be of the same length");
        }
        double result = 0.0;
        for(int i = 0; i < a.length; ++i) {
            result += a[i] * b[i];
        }
        return result;
    }

    public static double[] dot(double[][] a, double[] b) {
        double[] result = new double[a[0].length];
        for(int i = 0; i < a.length; ++i) {
            result[i] = dot(a[i], b);
        }
        return result;
    }

    // TODO check all b[0].length and a[0].length dimensions
    public static double[] dot(double[] a, double[][] b) {
        double[] result = new double[b[0].length];
        for(int i = 0; i < b.length; ++i) {
            if(a.length != b[i].length) {
                throw new IllegalArgumentException("Both arrays must be of the same length");
            }
            for(int j = 0; j < b[0].length; ++ j) {
                result[i] += a[j] * b[j][i];
            }
        }
        return result;
    }

    public static double[] dot(double[] a, double b) {
        return multiply(a, b);
    }

    public static void main(String[] args) {

        System.out.println(min(2,9,1,5));
        double[] a = {1, 2, 3};
        double[] b = {4, 5, 6};

        System.out.println(Arrays.toString(kron(a, b)));

        double[][] aa = {{1, 2}, {3, 4}};
        double[] bb = {5, 6};
        System.out.println(Arrays.toString(dot(bb, aa)));

        System.out.println(Arrays.toString(dot(aa, bb)));
    }
}
