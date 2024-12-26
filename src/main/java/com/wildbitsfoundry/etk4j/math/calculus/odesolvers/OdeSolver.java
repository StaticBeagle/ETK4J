package com.wildbitsfoundry.etk4j.math.calculus.odesolvers;

import com.wildbitsfoundry.etk4j.util.Tuples;

import java.util.Arrays;
// TODO document this class
// TODO implement 8
// TODO implement BDF
public abstract class OdeSolver {

    protected double t;
    protected Double tOld;
    protected double[] y;
    protected ODESystemOfEquations systemOfEquations;
    protected double tBound;
    protected double direction;
    protected String status; // TODO change to enum
    protected int n;

    public OdeSolver(ODESystemOfEquations systemOfEquations, double t0, double[] y0, Double tBound) {
        this.t = t0;
        if (!Arrays.stream(y0).allMatch(Double::isFinite)) {
            throw new IllegalArgumentException("All components of the initial state y0 must be finite.");
        }
        this.systemOfEquations = systemOfEquations;
        this.direction = (int) Math.signum(tBound - t0);
        this.y = y0;
        this.n = y0.length;
        this.status = "running";
        // nfev
        // njev
        // nlu
    }
// TODO test this
    private Double stepSize() {
        return tOld == null ? null : Math.abs(t - tOld);
    }

    public String step() {
        if (!this.status.equals("running")) {
            throw new RuntimeException("Attempt to step on a failed or finished solver.");
        }

        if (this.n == 0 || this.t == this.tBound) {
            this.tOld = this.t;
            this.t = this.tBound;
            this.status = "finished";
            return null;
        }
        Tuples.Tuple2<Boolean, String> result = stepImpl();

        if (!result.getItem1()) {
            this.status = "failed";
        } else {
            this.tOld = this.t;
            if (this.direction * (this.t - this.tBound) >= 0) {
                status = "finished";
            }
        }
        return result.getItem2();
    }

    protected abstract Tuples.Tuple2<Boolean, String> stepImpl();

    // TODO dense output
//    public Object getDenseOutput() {
//        if (tOld == null) {
//            throw new RuntimeException("Dense output is available after a successful step was made.");
//        }
//
//        if (n == 0 || t == tOld) {
//            // return ConstantDenseOutput(told, t, y);
//            return null; // TODO
//        } else {
//            return getDenseOutputImpl();
//        }
//    }
//
//    protected abstract DenseOutput getDenseOutputImpl();
//
//    public static class DenseOutput {
//
//    }
//
//    public static class ConstantDenseOutput extends DenseOutput {
//
//    }
}
