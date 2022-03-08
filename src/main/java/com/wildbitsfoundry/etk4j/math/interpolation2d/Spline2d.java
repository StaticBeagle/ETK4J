package com.wildbitsfoundry.etk4j.math.interpolation2d;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.functions.BivariateFunction;
import com.wildbitsfoundry.etk4j.math.interpolation.CubicSpline;
import com.wildbitsfoundry.etk4j.math.interpolation.Interpolation;
import com.wildbitsfoundry.etk4j.math.interpolation.LinearSpline;
import com.wildbitsfoundry.etk4j.math.interpolation.Spline;
// TODO document
/**
 *
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

    public static Spline2d newBicubicSpline(double[] x, double[] y, double[][] z) {
        final int rows = y.length;
        final int cols = x.length;
        final int order = 4;

        if (cols % order != 0) {
            throw new IllegalArgumentException(String.format("x length has to be a multiple of %d", order));
        }
        if (rows % order != 0) {
            throw new IllegalArgumentException(String.format("y length has to be a multiple of %d", order));
        }
        if (z.length != y.length) {
            throw new IllegalArgumentException(
                    String.format("The number of arrays in z has to be a multiple of %d", order));
        }
        for(int i = 0; i < z.length; ++i) {
            if(cols != z[i].length) {
                throw new IllegalArgumentException("The number of values in each array in z has to be the same");
            }
        }
        double[] yt = Arrays.copyOf(y, rows);
        Spline[] splines = new Spline[rows];
        for (int i = 0; i < rows; ++i) {
            if (cols != z[i].length) {
                throw new IllegalArgumentException("all arrays in z must have the same length");
            }
            splines[i] = CubicSpline.newCubicSplineInPlace(x, z[i]);
        }
        return new Spline2d(yt, splines, order);
    }

    public static Spline2d newBilinearSpline(double[] x, double[] y, double[][] z) {
        final int rows = y.length;
        final int cols = x.length;
        final int order = 2;

        if (cols % order != 0) {
            throw new IllegalArgumentException(String.format("x length has to be a multiple of %d", order));
        }
        if (rows % order != 0) {
            throw new IllegalArgumentException(String.format("y length has to be a multiple of %d", order));
        }
        if (z.length != y.length) {
            throw new IllegalArgumentException(
                    String.format("The number of arrays in z has to be a multiple of %d", order));
        }
        for(int i = 0; i < z.length; ++i) {
            if(cols != z[i].length) {
                throw new IllegalArgumentException("The number of values in each array in z has to be the same");
            }
        }
        double[] yt = Arrays.copyOf(y, rows);
        Spline[] splines = new Spline[rows];
        for (int i = 0; i < rows; ++i) {
            if (cols != z[i].length) {
                throw new IllegalArgumentException("all arrays in z must have the same length");
            }
            splines[i] = LinearSpline.newLinearSplineInPlace(x, z[i]);
        }
        return new Spline2d(yt, splines, order);
    }

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
