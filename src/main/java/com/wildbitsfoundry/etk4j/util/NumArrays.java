package com.wildbitsfoundry.etk4j.util;

import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.IntStream;

import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.interpolation.CubicSpline;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;

/**
 * The {@code NumArrays} class provides methods to manipulate arrays of native {@code double} values.
 */
public final class NumArrays {

    private NumArrays() {
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
        double[] result = Arrays.copyOf(a, a.length);
        subtractElementWiseInPlace(result, b);
        return result;
    }

    public static double[] subtract(double a, double[] b) {
        final int length = b.length;
        double[] result = new double[length];
        for (int i = 0; i < length; ++i) {
            result[i] = a - b[i];
        }
        return result;
    }

    public static double[] subtract(double[] a, double b) {
        final int length = a.length;
        double[] result = new double[length];
        for (int i = 0; i < length; ++i) {
            result[i] = a[i] - b;
        }
        return result;
    }

    public static void subtractElementWiseInPlace(double[] a, double[] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b dimensions must match.");
        }
        for (int i = 0; i < a.length; ++i) {
            a[i] = a[i] - b[i];
        }
    }

    public static double[] abs(double[] a) {
        double[] result = new double[a.length];
        for(int i = 0; i < a.length; ++i) {
            result[i] = Math.abs(a[i]);
        }
        return result;
    }

    public static double[] convolution(double[] a, double[] b) {
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

    // TODO these element-wise functions can be rewritten without the elementWise
    // I think the add function overlaps
    // Maybe intead of add in place, add equals?
    public static double[] addElementWise(double[] a, double b) {
        double[] result = Arrays.copyOf(a, a.length);
        addElementWiseInPlace(result, b);
        return result;
    }

    public static double[] addElementWise(double[] a, double[] b) {
        double[] result = Arrays.copyOf(a, a.length);
        addElementWiseInPlace(result, b);
        return result;
    }

    public static void addElementWiseInPlace(double[] a, double b) {
        for (int i = 0; i < a.length; ++i) {
            a[i] += b;
        }
    }

    public static void addElementWiseInPlace(double[] a, double[] b) {
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b dimensions must match");
        }
        for (int i = 0; i < a.length; ++i) {
            a[i] += b[i];
        }
    }

    public static double[] multiplyElementWise(double[] a, double d) {
        double[] result = Arrays.copyOf(a, a.length);
        multiplyElementWiseInPlace(result, d);
        return result;
    }

    public static void multiplyElementWiseInPlace(double[] a, double d) {
        final int length = a.length;
        for (int i = 0; i < length; ++i) {
            a[i] *= d;
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

    public static int[] argSort(double[] a) {
        Integer[] indexes = IntStream.range(0, a.length).boxed().toArray(Integer[]::new);
        Arrays.sort(indexes, Comparator.comparing(i -> a[i]));
        return Arrays.stream(indexes).mapToInt(Integer::intValue).toArray();
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
    public static double[] kronecker(final double[] a, final double[] b) {
        int aLength = a.length;
        int bLength = b.length;
        double[] c = new double[aLength * bLength];
        if (a.length == 1) {
            return NumArrays.multiplyElementWise(b, a[0]);
        }
        if (b.length == 1) {
            return NumArrays.multiplyElementWise(a, b[0]);
        }
        for (int i = 0; i < aLength; i++) {
            for (int j = 0; j < bLength; j++) {
                c[j + i * (bLength)] = a[i] * b[j];
            }
        }
        return c;
    }

    public static double[] concatenate(double[] a, double[] b) {
        int aLength = a.length;
        int bLength = b.length;
        double[] result = new double[aLength + bLength];
        System.arraycopy(a, 0, result, 0, aLength);
        System.arraycopy(b, 0, result, aLength, bLength);
        return result;
    }

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

    public static double[][] concatenate(double[][] a, double[][] b) {
        if(a.length != b.length) {
            throw new IllegalArgumentException("Both arrays must have the same number of rows.");
        }
        double[][] result = new double[a.length][];
        for(int i = 0; i < a.length; ++i) {
            result[i] = concatenate(a[i], b[i]);
        }
        return result;
    }

    public static double[][] stack(double[] a, double[] b) {
        double[][] result = new double[2][];
        result[0] = Arrays.copyOf(a, a.length);
        result[1] = Arrays.copyOf(b, b.length);
        return result;
    }

    public static double[][] stack(double[][] a, double[] b) {
        double[][] result = new double[a.length + 1][];
        for(int i = 0; i < a.length - 1; ++i) {
            result[i] = Arrays.copyOf(a[i], a[i].length);
        }
        result[result.length - 1] = Arrays.copyOf(b, b.length);
        return result;
    }

    public static double[][] stack(double[] a, double[][] b) {
        double[][] result = new double[b.length + 1][];
        result[0] = Arrays.copyOf(a, a.length);
        for(int i = 1; i < a.length; ++i) {
            result[i] = Arrays.copyOf(b[i], b[i].length);
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

    public static double sumSquares(double[] a) {
        final int length = a.length;
        double result = 0;
        for (int i = 0; i < length; ++i) {
            result += a[i] * a[i];
        }
        return result;
    }

    public static double[] cumulativeSum(double[] a) {
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

    public static double[] difference(double[] a) {
        double[] result = new double[a.length - 1];
        for(int i = 0; i < result.length; ++i) {
            result[i] = a[i + 1] - a[i];
        }
        return result;
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
        // TODO
        // check bounds
        // check that arrays have the same size

        double[] grad = new double[a.length];
        final int end = a.length - 1;
        grad[0] = (a[1] - a[0]) / (h[1] - h[0]); // once sided difference
        grad[end] = (a[end] - a[end - 1]) / (h[end] - h[end - 1]); // one sided difference

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

    public static double[] dot(double[] a, double[][] b) {
        if(a.length != b.length) {
            throw new IllegalArgumentException("The number of elements in a must match the number of rows in b.");
        }
        Matrix A = new Matrix(a, 1);
        A.multiplyEquals(new Matrix(b));
        return A.getArray();
    }

    public static double[] dot(double[] a, double b) {
        return multiplyElementWise(a, b);
    }

    public static double product(double[] a) {
        double prod = 1.0;
        for(int i = 0; i < a.length; ++i) {
            prod *= a[i];
        }
        return prod;
    }

    public static boolean allClose(double[] a, double target) {
        for(int i = 0; i < a.length; ++i) {
            if(!MathETK.isClose(a[i], target)) {
                return false;
            }
        }
        return true;
    }

    public static boolean allClose(double[] a, double target, double absTol, double relTol) {
        for(int i = 0; i < a.length; ++i) {
            if(!MathETK.isClose(a[i], target, absTol, relTol)) {
                return false;
            }
        }
        return true;
    }

    public static boolean allClose(double[] a, double target, double absTol) {
        for(int i = 0; i < a.length; ++i) {
            if(!MathETK.isClose(a[i], target, absTol)) {
                return false;
            }
        }
        return true;
    }

    public static double[] ones(int n) {
        double[] result = new double[n];
        Arrays.fill(result, 1.0);
        return result;
    }

    // TODO concat and stack operations in place?

    // https://stackoverflow.com/a/17634025/6383857
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

    public static void main(String[] args) {

        System.out.println(min(2,9,1,5));
        double[] a = {1, 2, 3};
        double[] b = {4, 5, 6};

        System.out.println(Arrays.toString(kronecker(a, b)));

        double[][] aa = {{1, 2}, {3, 4}};
        double[] bb = {5, 6};
        System.out.println(Arrays.toString(dot(bb, aa)));

        System.out.println(Arrays.toString(dot(aa, bb)));

        System.out.println(product(a));

        // TODO unit tests
        double[] aaa = {1, 1, 2, 3, 5, 8, 13, 21};
        System.out.println(Arrays.toString(difference(aaa)));

        double[] aaaa = {0, 1, 1, 0, 2, 3, 5, 0, 8, 13, 21};
        System.out.println(argMin(aaaa));

        System.out.println(Arrays.toString(add(a, b)));
        System.out.println(Arrays.toString(addElementWise(a, b)));

        System.out.println(Arrays.toString(convolution(aaa, aaaa)));

        System.out.println(Arrays.toString(convolution(aaa, bb)));
    }
}
