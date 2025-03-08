package com.wildbitsfoundry.etk4j.math.calculus.odesolver;

import com.wildbitsfoundry.etk4j.util.DoubleArrays;

import java.util.Arrays;

public class RungeKuttaDenseOutput extends DenseOutput {

    private final double h;
    private final int order;
    private final double[] yOld;
    private final double[][] Q;

    public RungeKuttaDenseOutput(double tOld, double t, double[] yOld, double[][] Q) {
        super(tOld, t);
        this.h = t - tOld;
        this.order = Q[0].length - 1;
        this.yOld = yOld;
        this.Q = Q;
    }

    @Override
    public double[][] evaluateAt(double[] t) {
        double[] x = DoubleArrays.divideElementWise(DoubleArrays.subtractElementWise(t, tOld), h);
        double[][] p = new double[this.order + 1][x.length];
        for(int i = 0; i <= this.order; i++) {
            p[i] = Arrays.copyOf(x, x.length);
        }
        for(int j = 0; j < x.length; j++) {
            for(int i = 1; i < p.length; i++) {
                p[i][j] *= p[i - 1][j];
            }
        }
        double[][] y = DoubleArrays.dot(Q, p);
        DoubleArrays.multiplyElementWiseInPlace(y, h);
        for(int i = 0; i < y.length; i++) {
            DoubleArrays.addElementWiseInPlace(y[i], yOld[i]);
        }
        return y;
    }

}
