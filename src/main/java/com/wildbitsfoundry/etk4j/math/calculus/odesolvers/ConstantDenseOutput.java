package com.wildbitsfoundry.etk4j.math.calculus.odesolvers;

import java.util.Arrays;

public class ConstantDenseOutput extends DenseOutput {

    private final double[] value;

    public ConstantDenseOutput(double tOld, double t, double[] value) {
        super(tOld, t);
        this.value = value;
    }

    @Override
    public double[][] evaluateAt(double[] t) {
        double[][] result = new double[value.length][t.length];
        for(int i = 0; i < value.length; i++) {
            double[] tmp = new double[t.length];
            Arrays.fill(tmp, value[i]);
            result[i] = tmp;
        }
        return result;
    }
}
