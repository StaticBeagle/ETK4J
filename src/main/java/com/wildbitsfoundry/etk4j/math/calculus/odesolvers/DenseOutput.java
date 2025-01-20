package com.wildbitsfoundry.etk4j.math.calculus.odesolvers;

public abstract class DenseOutput {

    protected final double tOld;
    protected final double t;
    protected final double tMin;
    protected final double tMax;

    public DenseOutput(double tOld, double t) {
        this.tOld = tOld;
        this.t = t;
        this.tMin = Math.min(t, tOld);
        this.tMax = Math.max(t, tOld);
    }
}
