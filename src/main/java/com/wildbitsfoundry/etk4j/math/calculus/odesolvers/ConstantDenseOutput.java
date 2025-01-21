package com.wildbitsfoundry.etk4j.math.calculus.odesolvers;

public class ConstantDenseOutput extends DenseOutput {

    private final double[] value;

    public ConstantDenseOutput(double tOld, double t, double[] value) {
        super(tOld, t);
        this.value = value;
    }

    @Override
    public double[][] evaluateAt(double[] t) {
        return null;
    }
}
