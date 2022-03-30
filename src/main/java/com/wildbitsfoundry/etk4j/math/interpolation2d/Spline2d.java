package com.wildbitsfoundry.etk4j.math.interpolation2d;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.functions.BivariateFunction;
import com.wildbitsfoundry.etk4j.math.interpolation.CubicSpline;
import com.wildbitsfoundry.etk4j.math.interpolation.Interpolation;
import com.wildbitsfoundry.etk4j.math.interpolation.LinearSpline;
import com.wildbitsfoundry.etk4j.math.interpolation.Spline;

/**
 * The {@code Spline2d} class represents bicubic and bilinear interpolation using cubic and linear splines respectively.
 */
public class Spline2d implements BivariateFunction {

    private double[] y;
    private Spline[] splines;

    private int order;

    protected Spline2d(double[] y, Spline[] splines, int order) {
        this.y = y;
        this.order = order;
        this.splines = splines;
    }

    /**
     * Construct a bicubic spline using not-a-knot conditions.
     *
     * @param x The x coordinates of the spline.
     * @param y The y coordinates of the spline.
     * @param z The values of the spline at each {@code (x, y)}. Notice that for each x and y combination there should
     *          be a value in {@code z}. For example, consider: <br/>
     *
     * {@code x = [1, 2, 3, 4]}, {@code y = [1, 2, 3, 4]}, <br/>and {@code z(x, y) = x<sup>2</sup> * y<sup>2</sup>}. <br/>
     *          This means that the first row of {@code z} will be all the values obtained by fixing {@code x = 1} and
     *          iterating over all the values of {@code y}. <br />
     *          {@code z[0] = [z(1, 1), z(1, 2), z(1, 3), z(1, 4)}; <br />
     *          Similarly: <br />
     *          {@code z[1] = [z(2, 1), z(2, 2), z(2, 3), z(2, 4)}; <br />
     *          and so on. <br />
     *          Even though internally we iterate over {@code x} and {@code y} multiple times, only a single copy of
     *          each {@code (x, y)} is required while. Please refer to the example in {@link examples.Spline2dExample}.
     *
     * @return A bicubic spline.
     */
    public static Spline2d newBicubicSpline(double[] x, double[] y, double[][] z) {
        final int rows = y.length;
        final int cols = x.length;
        final int order = 4;

        if (cols % order != 0) {
            throw new IllegalArgumentException(String.format("x length has to be a multiple of %d.", order));
        }
        if (rows % order != 0) {
            throw new IllegalArgumentException(String.format("y length has to be a multiple of %d.", order));
        }
        if (z.length != y.length) {
            throw new IllegalArgumentException(
                    String.format("The number of arrays in z has to be a multiple of %d.", order));
        }

        double[] yt = Arrays.copyOf(y, rows);
        Spline[] splines = new Spline[rows];
        for (int i = 0; i < rows; ++i) {
            if (cols != z[i].length) {
                throw new IllegalArgumentException("The length of each array in z has to be the same.");
            }
            splines[i] = CubicSpline.newCubicSplineInPlace(x, z[i]);
        }
        return new Spline2d(yt, splines, order);
    }

    /**
     * Construct a bilinear spline using not-a-knot conditions.
     *
     * @param x The x coordinates of the spline.
     * @param y The y coordinates of the spline.
     * @param z The values of the spline at each {@code (x, y)}. Notice that for each x and y combination there should
     *          be a value in {@code z}. For example, consider: <br/>
     *
     * {@code x = [1, 2, 3, 4]}, {@code y = [1, 2, 3, 4]}, <br/>and {@code z(x, y) = x<sup>2</sup> * y<sup>2</sup>}. <br/>
     *          This means that the first row of {@code z} will be all the values obtained by fixing {@code x = 1} and
     *          iterating over all the values of {@code y}. <br />
     *          {@code z[0] = [z(1, 1), z(1, 2), z(1, 3), z(1, 4)}; <br />
     *          Similarly: <br />
     *          {@code z[1] = [z(2, 1), z(2, 2), z(2, 3), z(2, 4)}; <br />
     *          and so on. <br />
     *          Even though internally we iterate over {@code x} and {@code y} multiple times, only a single copy of
     *          each {@code (x, y)} is required while. Please refer to the example in {@link examples.Spline2dExample}.
     *
     * @return A bilinear spline.
     */
    public static Spline2d newBilinearSpline(double[] x, double[] y, double[][] z) {
        final int rows = y.length;
        final int cols = x.length;
        final int order = 2;

        if (cols % order != 0) {
            throw new IllegalArgumentException(String.format("x length has to be a multiple of %d.", order));
        }
        if (rows % order != 0) {
            throw new IllegalArgumentException(String.format("y length has to be a multiple of %d.", order));
        }
        if (z.length != y.length) {
            throw new IllegalArgumentException(
                    String.format("The number of arrays in z has to be a multiple of %d.", order));
        }

        double[] yt = Arrays.copyOf(y, rows);
        Spline[] splines = new Spline[rows];
        for (int i = 0; i < rows; ++i) {
            if (cols != z[i].length) {
                throw new IllegalArgumentException("The length of each array in z has to be the same.");
            }
            splines[i] = LinearSpline.newLinearSplineInPlace(x, z[i]);
        }
        return new Spline2d(yt, splines, order);
    }

    /**
     * Evaluate the spline.
     *
     * @param x The x coordinate at which to evaluate the spline.
     * @param y The y coordinate at which to evaluate the spline.
     * @return The value of the spline at {@code (x, y}.
     */
    @Override
    public double evaluateAt(double x, double y) {
        int index = this.findLeftIndex(y);

        double[] tmp = new double[order];
        for (int i = 0; i < order; ++i) {
            tmp[i] = splines[i + index].evaluateAt(x);
        }
        return Interpolation.spline(Arrays.copyOfRange(this.y, index, index + order), tmp, y);
    }

    protected int findLeftIndex(double y) {
        int index = Arrays.binarySearch(this.y, y);
        if (index >= 0) {
            return order * (index / order);
        }
        index = -(index + 2);
        boolean edge = (index + 1) % order == 0;
        return edge ? (index + 1) - (order >> 1) : order * (index / order);
    }
}
