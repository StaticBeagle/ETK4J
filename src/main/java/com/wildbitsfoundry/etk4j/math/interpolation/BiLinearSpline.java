package com.wildbitsfoundry.etk4j.math.interpolation;

import examples.Spline2DExample;

import java.util.Arrays;

public class BiLinearSpline extends Spline2D {

    protected BiLinearSpline(double[] y, Spline[] splines, int order) {
        super(y, splines, order);
    }

    /**
     * Construct a bilinear spline using not-a-knot conditions.
     *
     * @param x The x coordinates of the spline.
     * @param y The y coordinates of the spline.
     * @param z The values of the spline at each {@code (x, y)}. Notice that for each x and y combination there should
     *          be a value in {@code z}. For example, consider: <br>
     *
     * {@code x = [1, 2, 3, 4]}, {@code y = [1, 2, 3, 4]}, <br>and {@code z(x, y) = x<sup>2</sup> * y<sup>2</sup>}. <br>
     *          This means that the first row of {@code z} will be all the values obtained by fixing {@code x = 1} and
     *          iterating over all the values of {@code y}. <br>
     *          {@code z[0] = [z(1, 1), z(1, 2), z(1, 3), z(1, 4)}; <br>
     *          Similarly: <br>
     *          {@code z[1] = [z(2, 1), z(2, 2), z(2, 3), z(2, 4)}; <br>
     *          and so on. <br>
     *          Even though internally we iterate over {@code x} and {@code y} multiple times, only a single copy of
     *          each {@code (x, y)} is required while. Please refer to the example in {@link Spline2DExample}.
     *
     * @return A bilinear spline.
     */
    public static Spline2D newBilinearSpline(double[] x, double[] y, double[][] z) {
        final int rows = y.length;
        final int cols = x.length;
        final int order = 2;

        if (z.length != y.length) {
            throw new IllegalArgumentException("The length of z has to be the same length as y.");
        }

        double[] yt = Arrays.copyOf(y, rows);
        Spline[] splines = new Spline[rows];
        for (int i = 0; i < rows; ++i) {
            if (cols != z[i].length) {
                throw new IllegalArgumentException("The length of each array in z has to be the same.");
            }
            splines[i] = LinearSpline.newLinearSplineInPlace(x, z[i]);
        }
        return new BiLinearSpline(yt, splines, order);
    }

}
