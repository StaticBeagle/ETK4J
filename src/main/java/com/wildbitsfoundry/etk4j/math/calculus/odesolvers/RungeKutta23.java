package com.wildbitsfoundry.etk4j.math.calculus.odesolvers;

import com.wildbitsfoundry.etk4j.math.functions.BivariateFunction;
// TODO document this class
public class RungeKutta23 extends RungeKutta {

    private static final double[][] A = {
            {0, 0, 0},
            {1 / 2d, 0, 0},
            {0, 3 / 4d, 0}
    };
    private static final double[] B = {2 / 9d, 1 / 3d, 4 / 9d};
    private static final double[] C = {0, 1 / 2d, 3 / 4d};
    private static final double[] E = {5 / 72d, -1 / 12d, -1 / 9d, 1 / 8d};
    private static final double[][] P = {
            {0, 1, -2 / 3d},
            {0, 4 / 3d, -8 / 9d},
            {0, -1, 1}
    };

    public RungeKutta23(ODESystemOfEquations systemOfEquations, double t0, double[] y0, Double tBound) {
        super(systemOfEquations, t0, y0, tBound, 2, 3, A, B, C, E, P);
    }

    public RungeKutta23(BivariateFunction func, double t0, double y0, Double tBound) {
        super(func, t0, y0, tBound, 2, 3, A, B, C, E, P);
    }

    public RungeKutta23(ODESystemOfEquations systemOfEquations, double t0, double[] y0, Double tBound, double maxStep,
                        double rTol, double aTol, Double firstStep) {
        super(systemOfEquations, t0, y0, tBound, maxStep, rTol, aTol, firstStep, 2, 3, A, B, C, E, P);
    }

    public RungeKutta23(BivariateFunction func, double t0, double y0, Double tBound, double maxStep,
                        double rTol, double aTol, Double firstStep) {
        super(func, t0, y0, tBound, maxStep, rTol, aTol, firstStep, 2, 3, A, B, C, E, P);
    }
}
