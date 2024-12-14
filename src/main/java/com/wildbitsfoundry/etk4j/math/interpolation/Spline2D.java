package com.wildbitsfoundry.etk4j.math.interpolation;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.functions.BivariateFunction;
import examples.Spline2DExample;

/**
 * The {@code Spline2D} class represents bicubic and bilinear interpolation using cubic and linear splines respectively.
 */
public abstract class Spline2D implements BivariateFunction {

    private double[] y;
    private Spline[] splines;

    private int order;

    protected Spline2D(double[] y, Spline[] splines, int order) {
        this.y = y;
        this.order = order;
        this.splines = splines;
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
        if(index == this.y.length - 1) {
            for (int i = this.y.length - order, j = 0; i < this.y.length; ++i, ++j) {
                tmp[j] = splines[i].evaluateAt(x);
            }
            return Interpolation.spline(Arrays.copyOfRange(this.y, this.y.length - order, order + 1), tmp, y);
        } else {
            for (int i = 0; i < order; ++i) {
                tmp[i] = splines[i + index].evaluateAt(x);
            }
            return Interpolation.spline(Arrays.copyOfRange(this.y, index, index + order), tmp, y);
        }
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
