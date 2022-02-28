package com.wildbitsfoundry.etk4j.math.optimize.solvers;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;

/**
 * The {@code Brent}'s class contains the Brent's root finding algorithm.
 * @see <a href="https://en.wikipedia.org/wiki/Brent%27s_method">Brent's method</a>
 */
public class Brent {

    private int maxNumberOfIterations = 100;
    private double absTol = 1e-9;
    private double relTol = 2 * ConstantsETK.DOUBLE_EPS;
    private double a;
    private double b;

    private UnivariateFunction function;


    /**
     * Constructs an instance of the {@code Brent}'s root finding algorithm.
     * @param function The function whose root is to be found.
     * @param a The lower end of the bracketing interval [a, b].
     * @param b The upper end of the bracketing interval [a, b].
     */
    public Brent(UnivariateFunction function, double a, double b) {
        this.function = function;
        this.a = a;
        this.b = b;
        if (b <= a) {
            throw new IllegalArgumentException("The upper bracket (b) must be greater than the lower bracket (a).");
        }
    }

    /**
     * Maximum number of iterations.
     * @param limit The maximum number of iterations allowed.
     */
    public Brent iterationLimit(int limit) {
        maxNumberOfIterations = limit;
        return this;
    }

    /**
     * Absolute tolerance.
     * @param tol The maximum allowed absolute tolerance.
     */
    public Brent absTolerance(double tol) {
        absTol = tol;
        return this;
    }

    /**
     * Relative tolerance.
     * @param tol The maximum allowed relative tolerance. This value must be bigger than 4 * {@link ConstantsETK#DOUBLE_EPS}.
     * @return
     */
    public Brent relTolerance(double tol) {
        if (tol < 4.0 * ConstantsETK.DOUBLE_EPS) {
            throw new IllegalArgumentException("Relative tolerance cannot be smaller than 4 * ConstantsETK.DOUBLE_EPS");
        }
        relTol = tol;
        return this;
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
     */

    /**
     * Find the root.
     * @return The {@link SolverResults} containing the root and other solver results.
     */
    SolverResults<Double> solve() {
        double xpre = a, xcur = b;
        double xblk = 0., fpre, fcur, fblk = 0., spre = 0., scur = 0., sbis;
        /* the tolerance is 2*delta */
        double delta;
        double stry, dpre, dblk;
        int currentIteration = 0;

        fpre = function.evaluateAt(xpre);
        fcur = function.evaluateAt(xcur);
        if (fpre * fcur > 0) {
            double error = Math.abs(Double.NaN);
            SolverResults<Double> solverResults = new SolverResults<>();
            solverResults.setSolverStatus("Root is not bracketed.");
            solverResults.setHasConverged(false);
            solverResults.setError(error);
            solverResults.setValue(xcur);
            solverResults.setNumberOfIterations(currentIteration);
            return solverResults;
        }
        if (fpre == 0) {
            double error = 0.0;
            SolverResults<Double> solverResults = new SolverResults<>();
            solverResults.setSolverStatus("Converged");
            solverResults.setHasConverged(true);
            solverResults.setError(error);
            solverResults.setValue(xpre);
            solverResults.setNumberOfIterations(currentIteration);
            return solverResults;
        }
        if (fcur == 0) {
            double error = 0.0;
            SolverResults<Double> solverResults = new SolverResults<>();
            solverResults.setSolverStatus("Converged");
            solverResults.setHasConverged(true);
            solverResults.setError(error);
            solverResults.setValue(xcur);
            solverResults.setNumberOfIterations(currentIteration);
            return solverResults;
        }

        //solver_stats->iterations = 0;
        for (currentIteration = 1; currentIteration <= maxNumberOfIterations; currentIteration++) {
            //solver_stats->iterations++;
            if (fpre != 0 && fcur != 0 &&
                    (Math.signum(fpre) != Math.signum(fcur))) {
                xblk = xpre;
                fblk = fpre;
                spre = scur = xcur - xpre;
            }
            if (Math.abs(fblk) < Math.abs(fcur)) {
                xpre = xcur;
                xcur = xblk;
                xblk = xpre;

                fpre = fcur;
                fcur = fblk;
                fblk = fpre;
            }

            delta = (absTol + relTol * Math.abs(xcur)) / 2;
            sbis = (xblk - xcur) / 2;
            if (fcur == 0 || Math.abs(sbis) < delta) {
                double error = Math.abs(xcur - xpre);
                SolverResults<Double> solverResults = new SolverResults<>();
                solverResults.setSolverStatus("Converged");
                solverResults.setHasConverged(true);
                solverResults.setError(error);
                solverResults.setValue(xcur);
                solverResults.setNumberOfIterations(currentIteration);
                return solverResults;
            }

            if (Math.abs(spre) > delta && Math.abs(fcur) < Math.abs(fpre)) {
                if (xpre == xblk) {
                    /* interpolate */
                    stry = -fcur * (xcur - xpre) / (fcur - fpre);
                } else {
                    /* extrapolate */
                    dpre = (fpre - fcur) / (xpre - xcur);
                    dblk = (fblk - fcur) / (xblk - xcur);
                    stry = -fcur * (fblk * dblk - fpre * dpre)
                            / (dblk * dpre * (fblk - fpre));
                }
                if (2 * Math.abs(stry) < Math.min(Math.abs(spre), 3 * Math.abs(sbis) - delta)) {
                    /* good short step */
                    spre = scur;
                    scur = stry;
                } else {
                    /* bisect */
                    spre = sbis;
                    scur = sbis;
                }
            } else {
                /* bisect */
                spre = sbis;
                scur = sbis;
            }

            xpre = xcur;
            fpre = fcur;
            if (Math.abs(scur) > delta) {
                xcur += scur;
            } else {
                xcur += (sbis > 0 ? delta : -delta);
            }

            fcur = function.evaluateAt(xcur);
        }
        double error = Math.abs(xcur - xpre);
        SolverResults<Double> solverResults = new SolverResults<>();
        solverResults.setSolverStatus("Maximum number of iterations exceeded");
        solverResults.setHasConverged(false);
        solverResults.setError(error);
        solverResults.setValue(xcur);
        solverResults.setNumberOfIterations(currentIteration);
        return solverResults;
    }
}
