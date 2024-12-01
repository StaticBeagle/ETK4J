package com.wildbitsfoundry.etk4j.math.calculus.odesolvers;

import com.wildbitsfoundry.etk4j.math.functions.BivariateFunction;

public class RungeKutta45 extends RungeKutta {

    private static final double[][] A = {
            {0, 0, 0, 0, 0},
            {1 / 5d, 0, 0, 0, 0},
            {3 / 40d, 9 / 40d, 0, 0, 0},
            {44 / 45d, -56 / 15d, 32 / 9d, 0, 0},
            {19372 / 6561d, -25360 / 2187d, 64448 / 6561d, -212 / 729d, 0},
            {9017 / 3168d, -355 / 33d, 46732 / 5247d, 49 / 176d, -5103 / 18656d}
    };
    private static final double[] B = {35 / 384d, 0, 500 / 1113d, 125 / 192d, -2187 / 6784d, 11 / 84d};
    private static final double[] C = {0, 1 / 5d, 3 / 10d, 4 / 5d, 8 / 9d, 1};
    private static final double[] E = {-71 / 57600d, 0, 71 / 16695d, -71 / 1920d, 17253 / 339200d, -22 / 525d, 1 / 40d};
    private static final double[][] P = {
            {1, -8048581381d / 2820520608d, 8663915743d / 2820520608d, -12715105075d / 11282082432d},
            {0, 0, 0, 0},
            {0, 131558114200d / 32700410799d, -68118460800d / 10900136933d, 87487479700d / 32700410799d},
            {0, -1754552775d / 470086768d, 14199869525d / 1410260304d, -10690763975d / 1880347072d},
            {0, 127303824393d / 49829197408d, -318862633887d / 49829197408d, 701980252875d / 199316789632d},
            {0, -282668133 / 205662961d, 2019193451 / 616988883d, -1453857185 / 822651844d},
            {0, 40617522 / 29380423d, -110615467 / 29380423d, 69997945 / 29380423d}
    };

    public RungeKutta45(ODESystemOfEquations systemOfEquations, double t0, double[] y0, Double tBound, double maxStep,
                        double rTol, double aTol, Double firstStep) {
        super(systemOfEquations, t0, y0, tBound, maxStep, rTol, aTol, firstStep, 4, 6, A, B, C, E, P);
    }

    public RungeKutta45(BivariateFunction func, double t0, double y0, Double tBound, double maxStep,
                        double rTol, double aTol, Double firstStep) {
        super(func, t0, y0, tBound, maxStep, rTol, aTol, firstStep, 4, 6, A, B, C, E, P);
    }
}
