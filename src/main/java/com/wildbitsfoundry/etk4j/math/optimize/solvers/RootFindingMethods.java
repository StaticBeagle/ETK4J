package com.wildbitsfoundry.etk4j.math.optimize.solvers;

import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;

public final class RootFindingMethods {

    private RootFindingMethods() {
    }

    public static class Root {
        private boolean converged;
        private double value;
        private String message;
        private int iterations;

        private Root(boolean converged, double value, String message, int iterations) {
            this.converged = converged;
            this.value = value;
            this.message = message;
            this.iterations = iterations;
        }

        public boolean isConverged() {
            return converged;
        }

        public double getValue() {
            return value;
        }

        public String getMessage() {
            return message;
        }

        public int getIterations() {
            return iterations;
        }
    }

    /***
     *
     * @param function
     * @param a
     * @param b
     * @param tol
     * @param maxIter
     * @return
     */
    public static Root ridders(UnivariateFunction function, double a, double b, double tol, int maxIter) {
        double fa, fb, fc, fx, c, s, dx, x = 0.0, xold = 0.0;
        fa = function.evaluateAt(a);
        fb = function.evaluateAt(b);


        if (fa * fb > 0.0) {
            return new Root(false, Double.NaN, "Root is not bracketed", 0);
        }
        for (int i = 1; i <= maxIter; i++) {
            c = 0.5 * (a + b);
            fc = function.evaluateAt(c);
            s = Math.sqrt(fc * fc - fa * fb);
            if (s == 0.0) {
                return new Root(true, x, "Root found", i);
            }
            dx = (c - a) * fc / s;
            if (fa - fb < 0.0) {
                dx = -dx;
            }
            x = c + dx;
            fx = function.evaluateAt(x);
            if (i > 1) // Check if we are done
            {
                if (Math.abs(x - xold) < tol) {
                    return new Root(true, x, "Root found", i);
                }
            }
            xold = x; // Update the root
            // Re-bracket and continue
            if (fc * fx > 0.0) {
                if (fa * fx < 0.0) {
                    b = x;
                    fb = fx;
                } else {
                    a = x;
                    fa = fx;
                }
            } else {
                a = c;
                b = x;
                fa = fc;
                fb = fx;
            }
        }
        return new Root(false, x, "Max number of iterations exceeded", maxIter);
    }

    /***
     * Reference: https://rosettacode.org/wiki/Roots_of_a_function#Brent.27s_Method
     * @param function
     * @param x0
     * @param x1
     * @param tol
     * @param maxIter
     * @return
     */
    public static Root brents(UnivariateFunction function, double x0, double x1, double tol, int maxIter) {
        double a = x0;
        double b = x1;
        double fa = function.evaluateAt(a);    // calculated now to save function calls
        double fb = function.evaluateAt(b);    // calculated now to save function calls
        double fs = 0;        // initialize

        if (!(fa * fb < 0)) {
            return new Root(false, Double.NaN, "Root is not bracketed", 0);
        }

        if (Math.abs(fa) < Math.abs(b))    // if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
        {
            // swap a & b
            double tmp = a;
            a = b;
            b = tmp;

            // swap fa & fb
            tmp = fa;
            fa = fb;
            fb = tmp;
        }

        double c = a;            // c now equals the largest magnitude of the lower and upper bounds
        double fc = fa;            // precompute function evalutation for point c by assigning it the same value as fa
        boolean mflag = true;        // boolean flag used to evaluate if statement later on
        double s = 0;            // Our Root that will be returned
        double d = 0;            // Only used if mflag is unset (mflag == false)

        int iter;
        for (iter = 1; iter < maxIter; ++iter) {
            // stop if converged on root or error is less than tolerance
            if (Math.abs(b - a) < tol) {
                return new Root(true, s, "Method Converged", iter);
            }

            if (fa != fc && fb != fc) {
                // use inverse quadratic interopolation
                s = (a * fb * fc / ((fa - fb) * (fa - fc)))
                        + (b * fa * fc / ((fb - fa) * (fb - fc)))
                        + (c * fa * fb / ((fc - fa) * (fc - fb)));
            } else {
                // secant method
                s = b - fb * (b - a) / (fb - fa);
            }

            // checks to see whether we can use the faster converging quadratic && secant methods or if we need to use bisection
            if (((s < (3 * a + b) * 0.25) || (s > b)) ||
                    (mflag && (Math.abs(s - b) >= (Math.abs(b - c) * 0.5))) ||
                    (!mflag && (Math.abs(s - b) >= (Math.abs(c - d) * 0.5))) ||
                    (mflag && (Math.abs(b - c) < tol)) ||
                    (!mflag && (Math.abs(c - d) < tol))) {
                // bisection method
                s = (a + b) * 0.5;
                mflag = true;
            } else {
                mflag = false;
            }

            fs = function.evaluateAt(s);    // calculate fs
            d = c;        // first time d is being used (was not used on first iteration because mflag was set)
            c = b;        // set c equal to upper bound
            fc = fb;    // set f(c) = f(b)

            if (fa * fs < 0)    // fa and fs have opposite signs
            {
                b = s;
                fb = fs;    // set f(b) = f(s)
            } else {
                a = s;
                fa = fs;    // set f(a) = f(s)
            }

            if (Math.abs(fa) < Math.abs(fb)) // if magnitude of fa is less than magnitude of fb
            {
                // swap a and b
                double tmp = a;
                a = b;
                b = tmp;

                // make sure f(a) and f(b) are correct after swap
                tmp = fa;
                fa = fb;
                fb = tmp;
            }

        }
        return new Root(false, Double.NaN, "Max number of iterations exceeded", maxIter);
    }


//    double bisectionMethod(double a, double b, double tol, int maxIter)
//    {
    // TODO
//    }
}
