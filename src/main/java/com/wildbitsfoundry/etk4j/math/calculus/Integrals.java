package com.wildbitsfoundry.etk4j.math.calculus;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;

public final class Integrals {
    private Integrals() {
    }

    public static double[] cummulativeTrapz(double... a) {
        final int length = a.length;
        double[] result = new double[length];
        for (int i = 1; i < length; ++i) {
            result[i] = result[i - 1] + (a[i] + a[i - 1]) * 0.5;
        }
        return result;
    }

    /***
     * Unit spaced trapezoidal integration.
     * @param a array of values of the function.
     * @return the approximate integral of {@code a} with unit spacing.
     */
    public static double trapz(double... a) {
        final int length = a.length;
        double result = 0.0;
        for (int i = 1; i < length; ++i) {
            result += (a[i] + a[i - 1]) * 0.5;
        }
        return result;
    }
    
    /***
     * Evaluates the definite integral from a to b using the trapezoidal rule.
     * @param func the function to be integrated.
     * @param a the starting point of the integration.
     * @param b the end point of the integration.
     * @param n the number of steps.
     * @return the approximate definite integral of {@code func} from a to b. 
     */
    public static double trapz(UnivariateFunction func, double a, double b, int n) {
        if(b < a) {
        	throw new IllegalArgumentException("b must be greater than a.");
        }
        if(n <= 0) {
        	throw new IllegalArgumentException("n has to be greater than zero.");
        }
    	double h = (b - a) / n;
    	double sum = func.evaluateAt(a) + func.evaluateAt(b);
    	for(int i = 1; i < n; ++i) {
    		sum += 2.0 * func.evaluateAt(a + i * h);
    	}
    	return 0.5 * h * sum;
    }
    
    /***
     * Evaluates the definite integral from a to b using the Simpson's 1/3 rule.
     * @param func the function to be integrated.
     * @param a the starting point of the integration.
     * @param b the end point of the integration.
     * @param n the number of steps.
     * @return the approximate definite integral of {@code func} from a to b. 
     */
    public static double simpson(UnivariateFunction func, double a, double b, double n) {
        if (n % 2 != 0) {
            throw new IllegalArgumentException("The number of elements must be even");
        }

        double even = 0, odd = 0;
        double h = (b - a) / n;
        for (int i = 1; i < n; i = i + 2) {
            odd += func.evaluateAt(a + i * h);
        }
        for (int i = 2; i < n - 1; i = i + 2) {
            even += func.evaluateAt(a + i * h);
        }
        return h / 3.0 * (func.evaluateAt(a) + 2.0 * even + 4.0 * odd + func.evaluateAt(b));
    }

    public static double qadrat(UnivariateFunction func, double a, double b,
                                 double absTol, double relTol) {
        return qadrat(func, a, b, absTol, relTol, 50);
    }

    public static double qadrat(UnivariateFunction func, double a, double b,
                                 double absTol, double relTol, int maxEval) {
        double x, f0, f2, f3, f5, f6, f7, f9, f14, hmin, hmax, re, ae, result;

        RefInteger numEval = new RefInteger(0);

        hmax = (b - a) / 16.0;
        if (hmax == 0.0) return 0.0;
        re = relTol;
        ae = 2.0 * absTol / Math.abs(b - a);
        hmin = Math.abs(b - a) * re;
        x = a;
        f0 = func.evaluateAt(x);
        x = a + hmax;
        f2 = func.evaluateAt(x);
        x = a + 2.0 * hmax;
        f3 = func.evaluateAt(x);
        x = a + 4.0 * hmax;
        f5 = func.evaluateAt(x);
        x = a + 6.0 * hmax;
        f6 = func.evaluateAt(x);
        x = a + 8.0 * hmax;
        f7 = func.evaluateAt(x);
        x = b - 4.0 * hmax;
        f9 = func.evaluateAt(x);
        x = b;
        f14 = func.evaluateAt(x);
        result = lint(func, a, b, f0, f2, f3, f5, f6, f7, f9, f14,
                hmin, re, ae, numEval, maxEval) * 16.0;
        return result;
    }


    static private double lint(UnivariateFunction func,
                               double x0, double xn, double f0, double f2, double f3,
                               double f5, double f6, double f7, double f9, double f14,
                               double hmin, double re, double ae, RefInteger numEval, int maxEval) {
        /* this function is internally used by QADRAT */

        if (numEval.getValue() >= maxEval) {
            throw new RuntimeException("Maximum number of evaluations reached in GaussLobatto");
        }

        numEval.setValue(numEval.getValue() + 1);

        double x, v, w, h, xm, f1, f4, f8, f10, f11, f12, f13;

        xm = (x0 + xn) / 2.0;
        h = (xn - x0) / 32.0;
        x = xm + 4.0 * h;
        f8 = func.evaluateAt(x);
        x = xn - 4.0 * h;
        f11 = func.evaluateAt(x);
        x = xn - 2.0 * h;
        f12 = func.evaluateAt(x);
        v = 0.330580178199226 * f7 + 0.173485115707338 * (f6 + f8) +
                0.321105426559972 * (f5 + f9) + 0.135007708341042 * (f3 + f11) +
                0.165714514228223 * (f2 + f12) + 0.393971460638127e-1 * (f0 + f14);
        x = x0 + h;
        f1 = func.evaluateAt(x);
        x = xn - h;
        f13 = func.evaluateAt(x);
        w = 0.260652434656970 * f7 + 0.239063286684765 * (f6 + f8) +
                0.263062635477467 * (f5 + f9) + 0.218681931383057 * (f3 + f11) +
                0.275789764664284e-1 * (f2 + f12) + 0.105575010053846 * (f1 + f13) +
                0.157119426059518e-1 * (f0 + f14);
        if (Math.abs(v - w) < Math.abs(w) * re + ae || Math.abs(h) < hmin)
            return h * w;
        else {
            x = x0 + 6.0 * h;
            f4 = func.evaluateAt(x);
            x = xn - 6.0 * h;
            f10 = func.evaluateAt(x);
            v = 0.245673430093324 * f7 + 0.255786258286921 * (f6 + f8) +
                    0.228526063690406 * (f5 + f9) + 0.500557131525460e-1 * (f4 + f10) +
                    0.177946487736780 * (f3 + f11) + 0.584014599347449e-1 * (f2 + f12) +
                    0.874830942871331e-1 * (f1 + f13) +
                    0.189642078648079e-1 * (f0 + f14);
            return ((Math.abs(v - w) < Math.abs(v) * re + ae) ? h * v :
                    (lint(func, x0, xm, f0, f1, f2, f3, f4, f5, f6, f7, hmin, re, ae, numEval, maxEval) -
                            lint(func, xn, xm, f14, f13, f12, f11, f10, f9, f8, f7, hmin, re, ae, numEval, maxEval)));
        }
    }

    private static class RefInteger {
        private int val;

        public RefInteger(int val) {
            this.val = val;
        }

        public int getValue() {
            return val;
        }

        public void setValue(int val) {
            this.val = val;
        }
    }

    public static void main(String[] args) {
        UnivariateFunction fx = x -> Math.sin(x);
        double[] e = new double[3];
        e[0] = ConstantsETK.FLOAT_EPS; // relative tol
        e[1] = ConstantsETK.FLOAT_EPS; // absolute tol
        double gg = qadrat(fx, 0, Math.PI / 2.0,  e[1], e[0]);
        System.out.println(gg);
        
        gg = trapz(fx, 0, Math.PI / 2.0, 1000);
        System.out.println(gg);
        
        gg = simpson(fx, 0, Math.PI / 2.0, 1000);
        System.out.println(gg);
        
        gg = trapz(new double[] {1, 4, 9, 16, 25, 36});
        System.out.println(gg);
        
        double[] ggg = cummulativeTrapz(new double[] {1, 4, 9, 16, 25});
        System.out.println(Arrays.toString(ggg));
    }
}
