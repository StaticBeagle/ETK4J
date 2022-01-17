package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import com.wildbitsfoundry.etk4j.math.optimize.solvers.NewtonRaphsonComplex;
import com.wildbitsfoundry.etk4j.math.optimize.solvers.Secant;
import com.wildbitsfoundry.etk4j.math.optimize.solvers.NewtonSolverResults;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;
import com.wildbitsfoundry.etk4j.util.NumArrays;
import static com.wildbitsfoundry.etk4j.signals.filters.Filters.*;

import java.util.Arrays;
import java.util.function.Function;

public class Bessel extends AnalogFilter {
    public enum FrequencyNormalization {
        PHASE,
        DELAY,
        MAGNITUDE
    }

    public static ZeroPoleGain besselap(int n) {
        return besselap(n, FrequencyNormalization.PHASE);
    }

    /*
        Port from scipy
        ``phase``
            The filter is normalized such that the phase response reaches its
            midpoint at an angular (e.g., rad/s) cutoff frequency of 1. This
            happens for both low-pass and high-pass filters, so this is the
            "phase-matched" case. [6]_
            The magnitude response asymptotes are the same as a Butterworth
            filter of the same order with a cutoff of `Wn`.
            This is the default, and matches MATLAB's implementation.
        ``delay``
            The filter is normalized such that the group delay in the passband
            is 1 (e.g., 1 second). This is the "natural" type obtained by
            solving Bessel polynomials
        ``mag``
            The filter is normalized such that the gain magnitude is -3 dB at
            angular frequency 1. This is called "frequency normalization" by
            Bond. [1]_
     */
    public static ZeroPoleGain besselap(int n, FrequencyNormalization norm) { // change norm to enum
        if(n == 0) {
            return new ZeroPoleGain(new Complex[]{}, new Complex[]{}, 1.0);
        }
        Complex[] zeros = {};
        Complex[] poles = Arrays.stream(besselZeros(n)).map(Complex::invert).toArray(Complex[]::new);
        double k = 1.0;

        double aLast = Math.floor(fallingFactorial(2 * n, n) / Math.pow(2, n));
        switch (norm) {
            case DELAY:
                k = aLast;
                break;
            case MAGNITUDE:
                k = aLast;
                double normFactor = normFactor(poles, k);
                for(int i = 0; i < n; ++i) {
                    poles[i].divideEquals(normFactor);
                }
                k = Math.pow(normFactor, -n) * aLast;
                break;
            case PHASE:
                for(int i = 0; i < n; ++i) {
                    poles[i].multiplyEquals(Math.pow(10, -Math.log10(aLast) / n));
                }
                break;
        };
        return new ZeroPoleGain(zeros, poles, k);
    }

    // TODO add normalization phase,delay,mag as an overload
    public static TransferFunction newLowPass(int n, double wn) {
        ZeroPoleGain zpk = besselap(n);
        return lpTolp(zpk, wn);
    }

    public static TransferFunction newHighPass(int n, double wn) {
        ZeroPoleGain zpk = besselap(n);
        return lpTohp(zpk, wn);
    }

    // TODO create exceptions. add checks to other filters
    public static TransferFunction newBandPass(int n, double wp1, double wp2) {
        if(n <= 0) {
            // throw
        }
        if(wp1 <= 0 || wp2 <= 0) {
            // throw
        }
        if(wp1 <= wp2) {
            // throw
        }
        ZeroPoleGain zpk = besselap(n);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpTobp(zpk, w0, bw);
    }

    public static TransferFunction newBandStop(int n, double wp1, double wp2) {
        ZeroPoleGain zpk = besselap(n);
        double w0 = Math.sqrt(wp1 * wp2);
        double bw = wp2 - wp1;
        return lpTobs(zpk, w0, bw);
    }

    private static Complex[] besselZeros(int n) {
        if (n == 0) {
            return new Complex[0];
        }
        Complex[] x0 = camposZeros(n);
        Complex[] x = aberth(fx, fp, x0, 1e-15, 50);

        for (int i = 0; i < n; ++i) {
            NewtonSolverResults<Complex> sr = new NewtonRaphsonComplex(
                    z -> {
                        Complex[] result = new Complex[1];
                        com.wildbitsfoundry.etk4j.math.specialfunctions.
                                Bessel.nonexpbesska01(n + 0.5, z.invert(), result, new Complex[1]);
                        return result[0];
                    },
                    z -> {
                        Complex[] a0 = new Complex[1];
                        com.wildbitsfoundry.etk4j.math.specialfunctions.
                                Bessel.nonexpbesska01(n - 0.5, z.invert(), a0, new Complex[1]);
                        a0[0] = a0[0].divide(z.pow(2).multiply(2));

                        Complex[] a1 = new Complex[1];
                        com.wildbitsfoundry.etk4j.math.specialfunctions.
                                Bessel.nonexpbesska01(n + 0.5, z.invert(), a1, new Complex[1]);
                        a1[0] = a1[0].divide(z.pow(2));

                        Complex[] a2 = new Complex[1];
                        com.wildbitsfoundry.etk4j.math.specialfunctions.
                                Bessel.nonexpbesska01(n + 1.5, z.invert(), a2, new Complex[1]);
                        a2[0] = a2[0].divide(z.pow(2).multiply(2));

                        return a0[0].subtract(a1[0]).add(a2[0]);
                    }, x[i])
                    .absTolerance(1e-15)
                    .relTolerance(0.0)
                    .iterationLimit(50)
                    .solve();
            x[i] = sr.getValue();
        }

        Complex[][] mean = new Complex[2][n];
        mean[0] = ComplexArrays.deepCopy(x);
        for(int i = 0; i < n; ++i) {
            mean[1][n - i - 1] = x[i].conj();
        }

        for(int j = 0; j < n; ++j) {
            Complex[] temp = new Complex[mean.length];
            for(int i = 0; i < mean.length; ++i) {
                temp[i] = mean[i][j];
            }
            x[j] = ComplexArrays.mean(temp);
        }

        // zeros should sum to -1
        if(ComplexArrays.sum(x).add(1.0).abs() > 1e-15) {
            // throw exception TODO
        }
        return x;
    }

    private static Function<Complex[], Complex[]> fx = x -> {
        final int length = x.length;
        Complex[][] result = new Complex[length][1];
        for (int i = 0; i < length; ++i) {
            com.wildbitsfoundry.etk4j.math.specialfunctions.
                    Bessel.nonexpbesska01(length + 0.5, x[i].invert(), result[i], new Complex[1]);
        }
        return ComplexArrays.flatten(result);
    };

    private static Function<Complex[], Complex[]> fp = x -> {
        final int length = x.length;
        Complex[][] result = new Complex[length][1];
        for (int i = 0; i < length; ++i) {
            Complex[] a0 = new Complex[1];
            com.wildbitsfoundry.etk4j.math.specialfunctions.
                    Bessel.nonexpbesska01(length - 0.5, x[i].invert(), a0, new Complex[1]);
            a0[0] = a0[0].divide(x[i].pow(2).multiply(2));

            Complex[] a1 = new Complex[1];
            com.wildbitsfoundry.etk4j.math.specialfunctions.
                    Bessel.nonexpbesska01(length + 0.5, x[i].invert(), a1, new Complex[1]);
            a1[0] = a1[0].divide(x[i].pow(2));

            Complex[] a2 = new Complex[1];
            com.wildbitsfoundry.etk4j.math.specialfunctions.
                    Bessel.nonexpbesska01(length + 1.5, x[i].invert(), a2, new Complex[1]);
            a2[0] = a2[0].divide(x[i].pow(2).multiply(2));

            result[i] = new Complex[]{a0[0].subtract(a1[0]).add(a2[0])};
        }
        return ComplexArrays.flatten(result);
    };

    private static Complex[] camposZeros(int n) {
        if (n == 1) {
            return new Complex[]{Complex.fromReal(-1.0)};
        }

        double s = Polynomial.polyval(new double[]{1, -3, 0, 2, 0, 0}, n);
        double b3 = Polynomial.polyval(new double[]{-8, 16}, n) / s;
        double b2 = Polynomial.polyval(new double[]{12, -12, -24}, n) / s;
        double b1 = Polynomial.polyval(new double[]{-2, -12, 24, 8}, n) / s;
        double b0 = Polynomial.polyval(new double[]{-1, 5, 0, -6, 0}, n) / s;

        double r = Polynomial.polyval(new double[]{1, 2, 0, 0}, n);
        double a1 = Polynomial.polyval(new double[]{-6, -6}, n) / r;
        double a2 = 6 / r;

        double[] k = NumArrays.linSteps(1, n);
        double[] x = Polynomial.polyval(new double[]{a2, a1, 0}, k);
        double[] y = Polynomial.polyval(new double[]{b3, b2, b1, b0}, k);
        return ComplexArrays.zip(x, y);
    }

    private static Complex[] aberth(Function<Complex[], Complex[]> fx, Function<Complex[], Complex[]> fp,
                                    Complex[] x0, double tol, int maxIter) {
        final int n = x0.length;
        Complex[] x = ComplexArrays.deepCopy(x0);

        Complex[] beta = new Complex[n];
        int iter = 0;
        while (iter < maxIter) {
            Complex[] num = fx.apply(x);
            Complex[] den = fp.apply(x);
            Complex[] alpha = new Complex[n];
            for (int i = 0; i < n; ++i) {
                alpha[i] = num[i].divide(den[i]).uminus();
            }

            for (int i = 0; i < n; ++i) {
                beta[i] = new Complex();
                for (int j = i + 1; j < n; ++j) {
                    beta[i].addEquals(x[i].subtract(x[j]).invert());
                }
                for (int j = 0; j < i; ++j) {
                    beta[i].addEquals(x[i].subtract(x[j]).invert());
                }
            }

            for (int i = 0; i < n; ++i) {
                x[i].addEquals(alpha[i].divide(alpha[i].multiply(beta[i]).add(1.0)));
            }

            // TODO
            // if not all are finite throw exception

            boolean done = true;
            for (int i = 0; i < n; ++i) {
                done &= alpha[i].abs() <= tol;
            }
            if (done) {
                break;
            }
            if (++iter >= maxIter) {
                // throw failed to converge
            }
        }
        return x;
    }

    private static double fallingFactorial(int x, int n) {
        double val = 1.0;
        for(int i = x - n + 1; i < x + 1; ++i) {
            val *= i;
        }
        return val;
    }

    private static double normFactor(Complex[] poles, double k) {
        UnivariateFunction gw = w -> {
            Complex prod = Complex.fromReal(1.0);
            for(int i = 0; i < poles.length; ++i) {
                prod.multiplyEquals(Complex.fromImaginary(w).subtract(poles[i]));
            }
            return prod.invert().multiply(k).abs();
        };
        UnivariateFunction cutoff = w -> gw.evaluateAt(w) - 1.0 / Math.sqrt(2.0);
        // 1.0 / Math.sqrt(2.0) = -3 db which is equal to 10 ^ (-3.0 / 20.0) TODO change -3 to an arbitrary input
        double result = Secant.solve(cutoff, 1.5, 1.5 * (1 + 1e-4), 1.48e-8, 0.0, 50);
        return result;
    }
}
