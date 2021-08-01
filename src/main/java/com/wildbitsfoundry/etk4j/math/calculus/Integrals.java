package com.wildbitsfoundry.etk4j.math.calculus;

import java.util.function.BiFunction;

import com.wildbitsfoundry.etk4j.exceptions.MaximumNumberOfEvaluationsReached;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import com.wildbitsfoundry.etk4j.util.RefValue.RefInteger;

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
     * Computes the approximate definite integral from a to b using the trapezoidal rule.
     * @param func the function to be integrated.
     * @param a the lower limit of the integral.
     * @param b the upper limit of the integral.
     * @param n the number of partitions.
     * @param params optional parameters passed to {@code func}.
     * @return the approximate definite integral of {@code func} from a to b. 
     */
    public static double trapz(BiFunction<Double, Object[], Double> func, double a, double b, int n, Object ...params) {
        if(b < a) {
        	throw new IllegalArgumentException("b must be greater than a.");
        }
        if(n <= 0) {
        	throw new IllegalArgumentException("n has to be greater than zero.");
        }
    	double h = (b - a) / n;
    	double sum = func.apply(a, params) + func.apply(b, params);
    	for(int i = 1; i < n; ++i) {
    		sum += 2.0 * func.apply(a + i * h, params);
    	}
    	return 0.5 * h * sum;
    }
    
    /***
     * Computes the approximate definite integral from a to b using the trapezoidal rule.
     * @param func the function to be integrated.
     * @param a the lower limit of the integral.
     * @param b the upper limit of the integral.
     * @param n the number of partitions.
     * @return the approximate definite integral of {@code func} from a to b. 
     */
    public static double trapz(UnivariateFunction func, double a, double b, int n) {
        return trapz((x, o) -> func.evaluateAt(x), a, b, n);
    }
   
    /***
     * Computes the approximate definite integral from a to b using the Simpson's 1/3 rule.
     * @param func the function to be integrated.
     * @param a the lower limit of the integral.
     * @param b the upper limit of the integral.
     * @param n the number of partitions. Must be even.
     * @param params optional parameters passed to {@code func}.
     * @return the approximate definite integral of {@code func} from a to b. 
     */
    public static double simpson(BiFunction<Double, Object[], Double> func, double a, double b, int n, Object ...params) {
    	if(n % 2 != 0) {
    		throw new IllegalArgumentException("The number of partitions (n) must be even");
    	}
        double even = 0, odd = 0;
        double h = (b - a) / n;
        for (int i = 1; i < n; i = i + 2) {
            odd += func.apply(a + i * h, params);
        }
        for (int i = 2; i < n - 1; i = i + 2) {
            even += func.apply(a + i * h, params);
        }
        return h / 3.0 * (func.apply(a, params) + 2.0 * even + 4.0 * odd + func.apply(b, params));
    }
    
    /***
     * Computes the approximate definite integral from a to b using the Simpson's 1/3 rule.
     * @param func the function to be integrated.
     * @param a the lower limit of the integral.
     * @param b the upper limit of the integral.
     * @param n the number of partitions. Must be even.
     * @param params optional parameters passed to {@code func}.
     * @return the approximate definite integral of {@code func} from a to b. 
     */
    public static double simpson(UnivariateFunction func, double a, double b, int n) {
        return simpson((x, o) -> func.evaluateAt(x), a, b, n);
    }
    
    
    /***
     * Computes the approximate definite integral from a to b using 1d quadrature.
     * @param func the function to be integrated.
     * @param a the lower limit of the integral.
     * @param b the upper limit of the integral.
     * @param absTol the absolute tolerance.
     * @param relTol the relative tolerance.
     * @param params optional parameters passed to {@code func}. 
     * @return the approximate definite integral of {@code func} from a to b. 
     */
    public static double qadrat(BiFunction<Double, Object[], Double> func, double a, double b,
                                 double absTol, double relTol, int maxEval, Object ...params) {
        double x, f0, f2, f3, f5, f6, f7, f9, f14, hmin, hmax, re, ae, result;

        RefInteger numEval = new RefInteger(0);

        hmax = (b - a) / 16.0;
        if (hmax == 0.0) return 0.0;
        re = relTol;
        ae = 2.0 * absTol / Math.abs(b - a);
        hmin = Math.abs(b - a) * re;
        x = a;
        f0 = func.apply(x, params);
        x = a + hmax;
        f2 = func.apply(x, params);
        x = a + 2.0 * hmax;
        f3 = func.apply(x, params);
        x = a + 4.0 * hmax;
        f5 = func.apply(x, params);
        x = a + 6.0 * hmax;
        f6 = func.apply(x, params);
        x = a + 8.0 * hmax;
        f7 = func.apply(x, params);
        x = b - 4.0 * hmax;
        f9 = func.apply(x, params);
        x = b;
        f14 = func.apply(x, params);
        result = lint(func, a, b, f0, f2, f3, f5, f6, f7, f9, f14,
                hmin, re, ae, numEval, maxEval, params) * 16.0;
        return result;
    }

    /***
     * Computes the approximate definite integral from a to b using 1d quadrature.
     * @param func the function to be integrated.
     * @param a the lower limit of the integral.
     * @param b the upper limit of the integral.
     * @param absTol the absolute tolerance.
     * @param relTol the relative tolerance.
     * @param 
     * @return the approximate definite integral of {@code func} from a to b. 
     */
    public static double qadrat(UnivariateFunction func, double a, double b,
                                 double absTol, double relTol, int maxEval) {
        return qadrat((x, o) -> func.evaluateAt(x), a, b, absTol, relTol, maxEval);
    }


    static private double lint(BiFunction<Double, Object[], Double> func,
                               double x0, double xn, double f0, double f2, double f3,
                               double f5, double f6, double f7, double f9, double f14,
                               double hmin, double re, double ae, RefInteger numEval, int maxEval, Object ...params) {
        /* this function is internally used by QADRAT */

        if (numEval.getValue() >= maxEval) {
        	String error = String.format("Maximum number of evaluations reached in qadrat.%n"
        			+ "If increasing the number of evaluations doesn't help, try loosening the%n"
        			+ "relative and absolute tolerances.");
            throw new MaximumNumberOfEvaluationsReached(error);
        }

        numEval.setValue(numEval.getValue() + 1);

        double x, v, w, h, xm, f1, f4, f8, f10, f11, f12, f13;

        xm = (x0 + xn) / 2.0;
        h = (xn - x0) / 32.0;
        x = xm + 4.0 * h;
        f8 = func.apply(x, params);
        x = xn - 4.0 * h;
        f11 = func.apply(x, params);
        x = xn - 2.0 * h;
        f12 = func.apply(x, params);
        v = 0.330580178199226 * f7 + 0.173485115707338 * (f6 + f8) +
                0.321105426559972 * (f5 + f9) + 0.135007708341042 * (f3 + f11) +
                0.165714514228223 * (f2 + f12) + 0.393971460638127e-1 * (f0 + f14);
        x = x0 + h;
        f1 = func.apply(x, params);
        x = xn - h;
        f13 = func.apply(x, params);
        w = 0.260652434656970 * f7 + 0.239063286684765 * (f6 + f8) +
                0.263062635477467 * (f5 + f9) + 0.218681931383057 * (f3 + f11) +
                0.275789764664284e-1 * (f2 + f12) + 0.105575010053846 * (f1 + f13) +
                0.157119426059518e-1 * (f0 + f14);
        if (Math.abs(v - w) < Math.abs(w) * re + ae || Math.abs(h) < hmin)
            return h * w;
        else {
            x = x0 + 6.0 * h;
            f4 = func.apply(x, params);
            x = xn - 6.0 * h;
            f10 = func.apply(x, params);
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
}
