package com.wildbitsfoundry.etk4j.math.calculus;

import com.wildbitsfoundry.etk4j.math.functions.BivariateFunction;
import com.wildbitsfoundry.etk4j.util.Tuples;

public abstract class OdeSolver {

    protected double t;
    protected Double tOld;
    protected double y;
    protected BivariateFunction func;
    protected double tBound;
    protected double direction;
    protected String status; // TODO change to enum
    protected int n;

    public OdeSolver(BivariateFunction func, double t0, double y0, Double tBound) {
        this.t = t0;
        if(!Double.isFinite(y0)) {
            throw new IllegalArgumentException("y0 must be a real value");
        }
        this.func = func;
        this.direction = (int) Math.signum(tBound - t0);
        this.y = y0;
        this.n = 1;
        this.status = "running";
        // nfev
        // njev
        // nlu
    }

    private Double stepSize() {
        return tOld == null ? null : Math.abs(t - tOld);
    }

    public String step() {
        if(status != "running") {
            throw new RuntimeException("Attempt to step on a failed or finished solver.");
        }

        if(n == 0 || t == tBound) {
            tOld = t;
            t = tBound;
            status = "finished";
            return null;
        }
        double t = this.t;
        Tuples.Tuple2<Boolean, String> result = stepImpl();

        if(!result.getItem1()) {
            status = "failed";
        } else {
            tOld = t;
            if(direction * (t - tBound) >= 0) {
                status = "finished";
            }
        }
        return result.getItem2();
    }

    protected abstract Tuples.Tuple2<Boolean, String> stepImpl();

    public Object getDenseOutout() {
        if(tOld == null) {
            throw new RuntimeException("Dense output is available after a successful step was made.");
        }

        if(n == 0 || t == tOld) {
            // return ConstantDenseOutput(told, t, y);
            return null; // TODO
        } else {
            return getDenseOutputImpl();
        }
    }

    protected abstract Object getDenseOutputImpl();

    public static class DenseOutput {

    }

    public static class ConstantDenseOutput extends DenseOutput {

    }
}
