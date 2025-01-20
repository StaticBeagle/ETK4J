package com.wildbitsfoundry.etk4j.math.calculus.odesolvers;

import com.wildbitsfoundry.etk4j.util.DoubleArrays;

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

    public double evaluateAt(double t) {
        double x = (t - tOld) / h;
        double[] p = DoubleArrays.repeat(new double[] {x}, this.order + 1);
        p = DoubleArrays.cumulativeProduct(p);
        double[] y = DoubleArrays.dot(Q, p);
        DoubleArrays.multiplyElementWiseInPlace(y, h);
        DoubleArrays.addElementWiseInPlace(y, yOld);
        return y[0];
    }

    public double[] evaluateAt(double[] t) {
        double[] x = DoubleArrays.divideElementWise(DoubleArrays.subtractElementWise(t, tOld), h);
        double[][] p = new double[this.order + 1][1];
        double[] pCol0 = new double[this.order + 1];
        double[] pCol1 = new double[this.order + 1];
        for(int i = 0; i <= this.order; i++) {
            pCol0[i] = x[0];
            pCol1[i] = x[1];
        }
        pCol0 = DoubleArrays.cumulativeProduct(pCol0);
        pCol1 = DoubleArrays.cumulativeProduct(pCol1);
        for(int i = 0; i <= this.order; i++) {
            p[i] = new double[] {pCol0[i], pCol1[i]};
        }
        double[] y = DoubleArrays.dot(Q, p)[0];
        DoubleArrays.multiplyElementWiseInPlace(y, h);
        if(yOld.length == 1) {
            DoubleArrays.addElementWiseInPlace(y, yOld[0]);
        } else {
            DoubleArrays.addElementWiseInPlace(y, yOld); // TODO needs to check if this need to be padded
        }
        return y;
    }

}
