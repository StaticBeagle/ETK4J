package com.wildbitsfoundry.etk4j.control;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrices;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;
import com.wildbitsfoundry.etk4j.util.NumArrays;

import java.util.Arrays;

/**
 * The {@code Transfer Function} class represents a continuous time system in the frequency domain. <br>
 * The internal structure of the class is a {@link RationalFunction}.
 */
public class TransferFunction extends LinearTimeInvariantSystem {
    private RationalFunction rf;

    /***
     * Constructs a {@code TransferFunction} from the given zeros and poles.The gain of the system is calculated from
     * the given zeros and poles.
     * @param zeros The zeros of the transfer function.
     * @param poles The poles of the transfer function.
     */
    public TransferFunction(Complex[] zeros, Complex[] poles) {
        rf = new RationalFunction(zeros, poles);
    }

    /***
     * Constructs a {@code TransferFunction} from the given zeros, poles, and gain.
     * @param zeros The zeros of the transfer function.
     * @param poles The poles of the transfer function.
     * @param gain The gain of the transfer function.
     */
    public TransferFunction(Complex[] zeros, Complex[] poles, double gain) {
        rf = new RationalFunction(zeros, poles, gain);
    }

    /**
     * Copy constructor.
     * @param tf The {@link TransferFunction} to be copied.
     */
    public TransferFunction(TransferFunction tf) {
        rf = new RationalFunction(tf.rf);
    }

    /**
     * Constructs a {@code TransferFunction} from the given numerator {@link Polynomial} and denominator
     * {@link Polynomial}. This constructor does not normalize the resulting {@link TransferFunction} in other words,
     * it doesn't force the leading coefficient of the denominator {@link Polynomial} to be one.
     * @param numerator Numerator polynomial.
     * @param denominator Denominator polynomial
     */
    public TransferFunction(Polynomial numerator, Polynomial denominator) {
        rf = new RationalFunction(numerator, denominator);
    }

    /**
     * Constructs a {@code TransferFunction} from the given numerator and denominator. The numerator and denominator
     * values represent the coefficients of the numerator and denominator polynomials respectively. The coefficients are
     * assumed to be in descending order i.e. {@code [1, 2, 3]} represents {@code x<sup>2</sup> + 2x + 3}.
     * This constructor normalizes resulting {@link TransferFunction} to force the leading coefficient of the
     * denominator {@link Polynomial} to be one.
     * @param numerator The coefficients of the numerator polynomial in descending order.
     * @param denominator The coefficients of the denominator polynomial in descending order.
     */
    public TransferFunction(double[] numerator, double[] denominator) {
        rf = new RationalFunction(numerator, denominator);
    }

    /**
     * Constructs a {@code Transfer Function} from the given {@link RationalFunction}. This constructor does not
     * normalize the input {@link RationalFunction} to force the leading coefficient of the denominator to be one and
     * the function is taken as is.
     * @param rf The rational function to be used internally by the transfer function.
     */
    public TransferFunction(RationalFunction rf) {
        this.rf = new RationalFunction(rf);
    }

    /**
     * Constructs a {@code Transfer Function} from the given {@link ZeroPoleGain} representation.
     * @param zpk The zero, pole, gain representation of the system.
     */
    public TransferFunction(ZeroPoleGain zpk) {
        this(zpk.getZeros(), zpk.getPoles(), zpk.getGain());
    }

    /**
     * {@link ZeroPoleGain} representation of the system.
     * @return The transfer function object broken into its zeros, poles, and gain.
     */
    @Override
    public ZeroPoleGain toZeroPoleGain() {
        Complex[] zeros = rf.getZeros();
        Complex[] poles = rf.getPoles();
        double k = rf.getNumerator().getCoefficientAt(0) / rf.getDenominator().getCoefficientAt(0);
        return new ZeroPoleGain(zeros, poles, k);
    }

    /**
     * Zeros of the system.
     * @return An array containing the zeros of the system.
     */
    public Complex[] getZeros() {
        return rf.getZeros();
    }

    /**
     * Poles of the system.
     * @return An array containing the poles of the system.
     */
    public Complex[] getPoles() {
        return rf.getPoles();
    }

    /**
     * Numerator {@code Polynomial}.
     * @return The numerator polynomial.
     */
    public Polynomial getNumerator() {
        return rf.getNumerator();
    }

    /**
     * Numerator coefficients.
     * @return An array containing the coefficients of the numerator in descending order.This is equivalent to
     * calling {@link #getNumerator()} and subsequently calling {@link Polynomial#getCoefficients()}.
     * However, this method pads the resulting array with zeros if the degree of the denominator is greater than the degree
     * of the numerator e.g.
     * <pre>
     *     numerator: x + 1             -> [0, 1, 1];<br>
     *     denominator: x^2 + 2 * x + 1 -> [1, 2, 1];
     * </pre>
     */
    public double[] getNumeratorCoefficients() {
        Polynomial num = rf.getNumerator();
        int numOrder = num.degree();
        int denOrder = rf.getDenominator().degree();
        if(numOrder >= denOrder) {
            return num.getCoefficients();
        }
        double[] result = new double[denOrder + 1];
        int start = denOrder - numOrder;
        System.arraycopy(num.getCoefficients(), 0, result, start, numOrder + 1);
        return result;
    }

    /**
     * Denominator {@link Polynomial}.
     * @return The denominator polynomial.
     */
    public Polynomial getDenominator() {
        return rf.getDenominator();
    }

    /**
     * Denominator coefficients.
     * @return An array containing the coefficients of the denominator in descending order. This is equivalent to
     * calling {@link #getDenominator()} and subsequently calling {@link Polynomial#getCoefficients()}.
     * However, this method pads the resulting array with zeros if the degree of the numerator is greater than the degree
     * of the denominator e.g.
     * <pre>
     *     numerator: x^2 + 2 * x + 1 -> [1, 2, 1];<br>
     *     denominator: x + 1         -> [0, 1, 1];
     * </pre>
     */
    public double[] getDenominatorCoefficients() {
        Polynomial den = rf.getDenominator();
        int numOrder = rf.getNumerator().degree();
        int denOrder = den.degree();
        if(denOrder >= numOrder) {
            return den.getCoefficients();
        }
        double[] result = new double[numOrder + 1];
        int start = numOrder - denOrder;
        System.arraycopy(den.getCoefficients(), 0, result, start, denOrder + 1);
        return result;
    }

    /**
     * Evaluate the system at a given frequency.
     * @param w The frequency at which to evaluate the system.
     * @return The complex frequency response of the system.
     */
    public Complex evaluateAt(double w) {
        return rf.evaluateAt(0.0, w);
    }

    /**
     * Magnitude of the system.
     * @param w Argument at which to evaluate the function.
     * @return The absolute value of the complex response.
     */
    public double calculateMagnitudeAt(double w) {
        return this.evaluateAt(w).abs();
    }

    /**
     * Magnitude of the system.
     * @param w Argument at which to evaluate the function.
     * @return The absolute value of the complex response.
     */
    public double[] calculateMagnitudeAt(double[] w) {
        double[] magnitude = new double[w.length];
        for (int i = 0; i < w.length; ++i) {
            magnitude[i] = this.calculateMagnitudeAt(w[i]);
        }
        return magnitude;
    }

    /***
     * Calculate the system wrapped phase response.
     * @param w the frequencies where the phase needs to be calculated at.
     * @return The phase response of the system in rad / s.
     */
    public double calculatePhaseAt(double w) {
        return this.evaluateAt(w).arg();
    }

    /***
     * Calculate the system wrapped phase response.
     * @param w the frequencies where the phase needs to be calculated at.
     * @return The phase response of the system in rad / s.
     */
    public double[] calculatePhaseAt(double[] w) {
        double[] phase = new double[w.length];
        for (int i = 0; i < w.length; ++i) {
            phase[i] = this.calculatePhaseAt(w[i]);
        }
        return phase;
    }

    /***
     * Calculate the system wrapped phase response.
     * @param w the frequency where the phase needs to be calculated at.
     * @return The phase response of the system in degrees.
     */
    public double calculatePhaseInDegreesAt(double w) {
        return Math.toDegrees(this.evaluateAt(w).arg());
    }

    /***
     * Calculate the system wrapped phase response.
     * @param w the frequencies where the phase needs to be calculated at.
     * @return The phase response of the system in degrees.
     */
    public double[] calculatePhaseInDegreesAt(double[] w) {
        double[] phase = new double[w.length];
        for (int i = 0; i < w.length; ++i) {
            phase[i] = this.calculatePhaseInDegreesAt(w[i]);
        }
        return phase;
    }

    /***
     * Calculates the phase at of the system a given frequency. </br>
     * This operation uses the zeros and poles of the system to calculate the phase as:
     * <pre>
     * 	&Sigma; Phase(Zeros) - &Sigma; Phase(Poles)
     * </pre>
     * If you need to calculate the magnitude and phase or just the phase for an array of frequencies,</br>
     * it might be more efficient to use {@link #evaluateAt(double)}, and then <strong>unwrap</strong> the
     * phase information</br> using {@link #unwrapPhase(double[])}.
     * @param w
     * @return the phase of the system in degrees
     */
    public double calculateUnwrappedPhaseInDegreesAt(double w) {
        Complex[] zeros = this.getZeros();
        Complex[] poles = this.getPoles();

        double phase = calculatePhase(zeros, w);
        phase -= calculatePhase(poles, w);

        return Math.toDegrees(phase);
    }

    /***
     * Calculate the phase for a given set of roots.
     * @param roots An array of roots to calculate the phase.
     * @param w The frequency at which to calculate the phase.
     * @return The phase of the system at the prescribed frequency.
     */
    private static double calculatePhase(Complex[] roots, double w) {
        double phase = 0.0;
        for (int i = 0; i < roots.length; ++i) {
            if (roots[i].imag() == 0) {
                // Single real root
                double alpha = -roots[i].real();
                phase += Math.atan2(w, alpha);
            } else {
                // Complex pair of roots
                double alpha = -roots[i].real();
                double beta = roots[i].imag();
                double fn = MathETK.hypot(alpha, beta);
                double zeta = alpha / fn;
                double ratio = w / fn;
                phase += Math.atan2(2.0 * zeta * ratio, 1 - ratio * ratio);
                // Roots are in pairs. Skip the conjugate since it was
                // already taken into account in the calculation above
                ++i;
            }
        }
        return phase;
    }

    /**
     * Frequency response of the system. The default value for the number of points is 200.
     * @return The complex frequency response of the system and the frequency at which each point was calculated.
     */
    public FrequencyResponse calculateFrequencyResponse() {
        return this.calculateFrequencyResponse(200);
    }

    /**
     * Frequency response of the system.
     * @param numberOfPoints The number of points in which to evaluate the system.
     * @return The complex frequency response of the system and the frequency at which each point was calculated.
     */
    public FrequencyResponse calculateFrequencyResponse(int numberOfPoints) {
        double[] frequencies = findFrequencies(numberOfPoints);
        return this.calculateFrequencyResponse(frequencies);
    }

    /**
     * Frequency response of the system.
     * @param w The frequencies at which to evaluate the system.
     * @return The frequency response of the system at each given frequency.
     */
    public FrequencyResponse calculateFrequencyResponse(double[] w) {
        Complex[] response = new Complex[w.length];
        for (int i = 0; i < w.length; ++i) {
            response[i] = this.evaluateAt(w[i]);
        }
        return new FrequencyResponse(response, w);
    }

    /**
     * Bode response of the system.
     * @return The magnitude in dB and phase in degrees of the system.
     */
    public BodeResponse calculateBode() {
        return this.calculateBode(200);
    }

    /**
     * Bode response of the system.
     * @param numberOfPoints The number of points in which to evaluate the system.
     * @return The magnitude in dB and phase in degrees of the system.
     */
    public BodeResponse calculateBode(int numberOfPoints) {
        double[] frequencies = findFrequencies(numberOfPoints);
        return this.calculateBode(frequencies);
    }

    /**
     * Bode response of the system.
     * @param w The frequencies at which to evaluate the system.
     * @return The magnitude in dB and phase in degrees of the system.
     */
    public BodeResponse calculateBode(double[] w) {
        double[] magnitudeIndB = new double[w.length];
        double[] phaseInDegrees = new double[w.length];
        for (int i = 0; i < w.length; ++i) {
            magnitudeIndB[i] = 20 * Math.log10(this.calculateMagnitudeAt(w[i]));
            phaseInDegrees[i] = this.evaluateAt(w[i]).arg() * (180 / Math.PI);
        }
        return new BodeResponse(magnitudeIndB, phaseInDegrees, w);
    }

    /**
     * Helper method to find the "interesting frequencies" where the system is exhibiting change.
     * @param numberOfPoints The number of points to calculate the frequency at.
     * @return
     */
    private double[] findFrequencies(int numberOfPoints) {
        Complex[] ep = this.getPoles();
        Complex[] tz = this.getZeros();

        // TODO check if length of the poles == 0?
        Complex[] ez = ComplexArrays.concatenate(
                Arrays.stream(ep).filter(c -> c.imag() >= 0.0).toArray(Complex[]::new),
                Arrays.stream(tz).filter(c -> c.abs() < 1e5 && c.imag() >= 0.0).toArray(Complex[]::new)
        );
        int[] integ = Arrays.stream(ez).mapToInt(c -> c.abs() < 1e-8 ? 1 : 0).toArray();
        double[] argument = new double[ez.length];
        for(int i = 0; i < argument.length; ++i) {
            argument[i] = 3.0 * Math.abs(ez[i].real() + integ[i]) + 1.5 * ez[i].imag();
        }
        int hiFreq = roundFrequency(Math.log10(NumArrays.max(argument)) + 0.5);

        for(int i = 0; i < argument.length; ++i) {
            argument[i] = Math.abs(ez[i].add(integ[i]).real()) + 2 * ez[i].imag();
        }
        int loFreq = roundFrequency(Math.log10(0.1 * NumArrays.min(argument)) - 0.5);

        return NumArrays.logSpace(loFreq, hiFreq, numberOfPoints);
    }

    private static int roundFrequency(double d) {
        // get numbers after the decimal point
        double decimal = d - Math.floor(d);
        if(Math.abs(decimal) == 0.5) {
            return (int) MathETK.roundEven(d);
        } else{
            return (int) Math.round(d);
        }
    }

    /**
     * Unwraps frequency in degrees. For example, consider a system where the phase response goes from 0°
     * to 360°. If the phase is wrapped, it will only vary between 180° and -180° e.g. (Careful, ASCII art coming your way)
     * <pre>
     *  Expected:                     Actual:
     *      |                           |
     *      |                           |
     *      |                           |
     *      ---------------------       ---------------------
     *    0 |*******                    |******* *******
     *  180 |       *                   |       *
     *  360 |       *******             |
     *  </pre>
     * @param phase The phase to unwrap. This operation is done in plase and the phase array contains the result of
     *              unwrapping the frequency.
     */
    public static void unwrapPhase(double[] phase) {
        int length = phase.length;
        double[] dp = new double[length];
        double[] dps = new double[length];
        double[] C = new double[length];
        double[] cumulativeSum = new double[length];

        double cutoff = 180.0;
        int j;

        // incremental phase variation
        for (j = 0; j < length - 1; j++) {
            dp[j] = phase[j + 1] - phase[j];
        }
        // equivalent phase variation in [-pi, pi]
        for (j = 0; j < length - 1; j++) {
            dps[j] = (dp[j] + 180.0) - Math.floor((dp[j] + 180.0) / (2 * 180.0)) * (2 * 180.0) - 180.0;
        }
        // preserve variation sign for +pi vs. -pi
        for (j = 0; j < length - 1; j++) {
            if ((dps[j] == -180.0) && (dp[j] > 0)) {
                dps[j] = 180.0;
            }
        }
        // incremental phase correction
        for (j = 0; j < length - 1; j++) {
            C[j] = dps[j] - dp[j];
        }
        // Ignore correction when incremental variation is smaller than cutoff
        for (j = 0; j < length - 1; j++) {
            if (Math.abs(dp[j]) < cutoff) {
                C[j] = 0;
            }
        }
        // Find cumulative sum of deltas
        cumulativeSum[0] = C[0];
        for (j = 1; j < length - 1; j++) {
            cumulativeSum[j] = cumulativeSum[j - 1] + C[j];
        }
        // Integrate corrections and add to P to produce smoothed phase values
        for (j = 1; j < length; j++) {
            phase[j] += cumulativeSum[j - 1];
        }
    }

    /**
     * {@code TransferFunction} addition.
     * @param tf The transfer function to add.
     * @return The sum of the two transfer functions.
     */
    public TransferFunction add(final TransferFunction tf) {
        return new TransferFunction(rf.add(tf.rf));
    }

    /**
     * {@code TransferFunction} addition.
     * @param d The scalar to add.
     * @return The sum of the transfer function and the scalar {@code d}.
     */
    public TransferFunction add(double d) {
        return new TransferFunction(rf.add(d));
    }

    /**
     * {@code TransferFunction} subtraction.
     * @param tf The transfer function to subtract.
     * @return The subtraction of the two transfer functions.
     */
    public TransferFunction subtract(TransferFunction tf) {
        return new TransferFunction(rf.subtract(tf.rf));
    }

    /**
     * {@code TransferFunction} multiplication.
     * @param tf The transfer function to multiply.
     * @return The multiplication of the two transfer functions.
     */
    public TransferFunction multiply(TransferFunction tf) {
        return new TransferFunction(rf.multiply(tf.rf));
    }

    /**
     * {@code TransferFunction} multiplication.
     * @param d The scalar to multiply.
     * @return The multiplication of the transfer function and the scalar {@code d}.
     */
    public TransferFunction multiply(double d) {
        return new TransferFunction(rf.multiply(d));
    }

    @Override
    public String toString() {
        String num = rf.getNumerator().toString();
        String den = rf.getDenominator().toString();

        int numLength = num.length();
        int denLength = den.length();

        int divLength = Math.max(numLength, denLength) + 2;
        String divider = String.format("%0" + divLength + "d", 0).replace('0', '-');

        int padLength = (int) Math.floor((divLength - Math.min(numLength, denLength)) * 0.5);
        String padding = String.format("%0" + padLength + "d", 0).replace('0', ' ');

        String[] format = null;

        if (numLength > denLength) {
            format = new String[]{" ", padding};
        } else if (numLength < denLength) {
            format = new String[]{padding, " "};
        } else {
            format = new String[]{" ", " "};
        }
        String nl = System.getProperty("line.separator");
        String result = format[0] + num + nl + divider + nl + format[1] + den;
        return result.replace('x', 's');
    }

    /**
     * Complex polynomial magnitude.
     * @param polynomial The polynomial to evaluate at (jw) and then calculate the magnitude.
     * @return {@code |P(jw)|}
     */
    private static Polynomial getPolynomialMagnitude(Polynomial polynomial) {
        double[] num = polynomial.getCoefficients();
        Complex[] numjw = evaluateAtjw(num);
        return new Polynomial(norm(numjw));
    }

    /***
     * Calculates the squared magnitude element wise
     * @param a Array of complex numbers to calculate the magnitude.
     * @return [norm0, norm1, .... normn]
     */
    private static double[] norm(Complex[] a) {
        Complex[] mag = ComplexArrays.convolution(a, ComplexArrays.conj(a));
        double[] coefficients = new double[mag.length];
        for (int i = 0; i < mag.length; ++i) {
            coefficients[i] = mag[i].real();
        }
        return coefficients;
    }

    /**
     * Evaluate the polynomial at s = jw.
     * @param poly The polynomial coefficients.
     * @return The result of evaluating the frequency domain function at s = jw e.g.
     * <pre>
     *     H(s) = s^2 + s + 1;<br>
     *     H(jw) = -jw^2 + jw + 1;
     * </pre>
     */
    private static Complex[] evaluateAtjw(double[] poly) {
        final int length = poly.length;
        Complex[] result = new Complex[length];

        for (int i = length - 1, j = 0; i >= 0; --i, ++j) {
            poly[i] *= getComplexSign(j);
            if (j % 2 == 0) {
                result[i] = Complex.fromReal(poly[i]);
            } else {
                result[i] = Complex.fromImaginary(poly[i]);
            }
        }
        return result;
    }

    /**
     * Complex sign.
     * @param power Power of the complex number j.
     * @return The result of j^power.
     */
    private static int getComplexSign(int power) {
        int result = power - 4 * (power / 4);
        return result == 2 || result == 3 ? -1 : 1;
    }

    /**
     * Calculate all frequencies where the gain crosses 0 dB.
     * @return An array containing all the frequencies where the gain (in dB) crosses 0 dB in ascending order.
     * <pre>
     *        |Num(s)| <br>
     * H(s) = -------- = 1 <br>
     *        |Den(s)| <br><br>
     *
     * |Num(s)| - |Den(s)| = 0 <br><br>
     *
     * Find roots of the resultant polynomial.
     * </pre>
     */
    public double[] calculateAllGainCrossoverFrequencies() {
        return this.calculateAllGainCrossoverFrequencies(Math.sqrt(ConstantsETK.DOUBLE_EPS));
    }

    /**
     * Calculate all frequencies where the gain crosses 0 dB.
     * @param tol The tolerance used to determine the smallest acceptable frequency that can be considered.
     * @return An array containing all the frequencies where the gain (in dB) crosses 0 in ascending order.
     * <pre>
     *        |Num(s)| <br>
     * H(s) = -------- = 1 <br>
     *        |Den(s)| <br><br>
     *
     * |Num(s)| - |Den(s)| = 0 <br><br>
     *
     * Find roots of the resultant polynomial.
     * </pre>
     */
    public double[] calculateAllGainCrossoverFrequencies(double tol) {

        Polynomial magPolyNum = getPolynomialMagnitude(rf.getNumerator());
        Polynomial magPolyDen = getPolynomialMagnitude(rf.getDenominator());

        Complex[] solution = magPolyNum.subtract(magPolyDen).calculateRoots();
        double[] wgc = new double[solution.length];
        int j = 0;
        for (int i = 0; i < solution.length; i++) {
            double real = solution[i].real();
            if (solution[i].imag() == 0 && real >= tol) {
                wgc[j] = real;
                j++;
            }
        }
        wgc = Arrays.copyOfRange(wgc, 0, j);
        Arrays.sort(wgc);
        return wgc;
    }

    /**
     * Gain cross over frequency. The tolerance used to discern the minimum acceptable frequency to be considered is
     * equal to: {@code Math.sqrt(DOUBLE_EPSILON)}.
     * @return The first frequency at where the gain (in dB) crosses 0 dB.
     */
    public double calculateGainCrossoverFrequency() {
        return calculateGainCrossoverFrequency(Math.sqrt(ConstantsETK.DOUBLE_EPS));
    }

    /**
     * Gain cross over frequency.
     * @tol The tolerance used to determine the smallest acceptable frequency that can be considered.
     * @return The first frequency at where the gain (in dB) crosses 0 dB.
     */
    public double calculateGainCrossoverFrequency(double tol) {
        double[] freqs = this.calculateAllGainCrossoverFrequencies(tol);
        // Get the first one where the magnitude was decreasing
        for (double f : freqs) {
            double slope = this.calculateMagnitudeAt(f) - this.calculateMagnitudeAt(f + 1e-12);
            if (slope > 0.0) {
                return f;
            }
        }
        return Double.NaN;
    }

    /**
     * Phase crossover frequency.
     * @return The first frequency at which the phase crosses -180°.
     */
    public double calculatePhaseCrossoverFrequency() {
        return this.calculatePhaseCrossoverFrequency(Math.sqrt(ConstantsETK.DOUBLE_EPS));
    }

    /**
     * Calculate all frequencies where the phase crosses 0 -180°.
     * @param tol The tolerance used to determine the smallest acceptable frequency that can be considered.
     * @return An array containing all the frequencies where the phase (in degrees) crosses -180° in ascending order.
     */
    public double calculatePhaseCrossoverFrequency(double tol) {
        double[] frequencies = this.calculateAllPhaseCrossoverFrequencies(tol);
        // Get the first one where the magnitude was decreasing
        for (double freq : frequencies) {
            double slope = this.calculateMagnitudeAt(freq) - this.calculateMagnitudeAt(freq + 1e-12);
            if (slope > 0.0) {
                return freq;
            }
        }
        return Double.NaN;
    }

    public double[] calculateAllPhaseCrossoverFrequencies() {
        return this.calculateAllPhaseCrossoverFrequencies(Math.sqrt(ConstantsETK.DOUBLE_EPS));
    }

    /**
     * Finds all the phase crossover frequencies by solving H(s) = H(-s). The phase crossover frequencies are defined
     * as the frequencies at which the phase of the system crosses -180°. Typically, the first frequency at which
     * the phase crosses -180° is the frequency of interest and in that case
     * {@link #calculatePhaseCrossoverFrequency()} can be used. <br>
     *
     * <pre>
     * Let's define our system as:
     *
     *         Num(jw)
     * H(jw) = -------
     *         Den(jw)
     *
     * Multiply by the conjugate of the denominator Den(jw)':
     *
     *         Num(jw)   Den(jw)'
     * H(jw) = ------- * -------
     *         Den(jw)   Den(jw)'
     *
     *  Now we equate the imaginary part of the numerator to zero and solve for the frequencies.
     *
     *   Nump(jw) = Num(jw) * Den(jw)'
     *   Nump(jw) = Real(Nump(jw)) + j * Imag(Nump(jw))
     *
     *   Solve for the roots of:
     *
     *   Real(Nump(jw)) = 0
     * </pre>
     * The resulting roots are all the phase crossover frequencies.
     * @return An array of the phase crossover frequencies in ascending order.
     */
    public double[] calculateAllPhaseCrossoverFrequencies(double tol) {
        Complex[] num = evaluateAtjw(rf.getNumerator().getCoefficients());
        Complex[] den = evaluateAtjw(rf.getDenominator().getCoefficients());

        double[] convLHS = NumArrays.convolution(ComplexArrays.imag(num), ComplexArrays.real(den));
        double[] convRHS = NumArrays.convolution(ComplexArrays.real(num), ComplexArrays.imag(den));

        Polynomial a = new Polynomial(convLHS);
        Polynomial b = new Polynomial(convRHS);
        a.subtractEquals(b);

        Complex[] roots = a.calculateRoots();
        double[] wpc = new double[roots.length];
        int j = 0;
        for (int i = 0; i < roots.length; i++) {
            double real = roots[i].real();
            if (roots[i].imag() == 0 && real >= tol) {
                wpc[j] = real;
                j++;
            }
        }
        wpc = Arrays.copyOfRange(wpc, 0, j);
        Arrays.sort(wpc);
        return wpc;
    }

    /**
     * Margins of the system.
     * @return
     * <ol>
     *     <li>Gain Margin</li>
     *     <li>Phase Margin </li>
     *     <li>Phase Crossover frequency.</li>
     *     <li>Gain Crossover frequency</li>
     * </ol>
     */
    public Margins calculateMargins() {
        double wcg = this.calculateGainCrossoverFrequency();
        double wcp = this.calculatePhaseCrossoverFrequency();

        double pm = Double.POSITIVE_INFINITY;
        if (!Double.isNaN(wcg)) {
            double phase = this.calculateUnwrappedPhaseInDegreesAt(wcg);
            pm = 180.0 + phase;
        }
        double gm = Double.POSITIVE_INFINITY;
        if (!Double.isNaN(wcp)) {
            double mag = this.calculateMagnitudeAt(wcp);
            gm = 1.0 / mag;
        }

        Margins m = new Margins(gm, pm, wcp, wcg);
        return m;
    }

    /** TODO
     * Substitute
     * @param d
     */
    public void substituteInPlace(double d) {
        rf.substituteInPlace(d);
    }

    /**
     * Proper {@TransferFunction}.
     * @return True if the order of the numerator is less or equal than the order of the denominator.
     */
    public boolean isProper() {
        return rf.isProper();
    }

    /**
     * Strictly proper {@code TransferFunction}.
     * @return True if the order of the numerator is less than the order of the denominator.
     */
    public boolean isStrictlyProper() {
        return rf.isStrictlyProper();
    }

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
     */
    /***
     * Transform SISO only single input single output TFs
     * @return
     * @throws ImproperTransferFunctionException If the order of the numerator is greater than the order of the
     * denominator.
     */
    @Override
    public StateSpace toStateSpace() {
        TransferFunction tf = new TransferFunction(this);

        tf.normalize();
        if (!tf.isProper()) {
            throw new ImproperTransferFunctionException("Transfer function must be proper.");
        }
        double[] num = tf.rf.getNumerator().getCoefficients();
        double[] den = tf.rf.getDenominator().getCoefficients();

        if (num.length == 0.0 || den.length == 0.0) {
            // Null system
            return new StateSpace();
        }
        // Pad numerator with zeros to match denominator size
        double[] numPadded = new double[den.length];
        System.arraycopy(num, 0, numPadded, den.length - num.length, num.length);

        double[][] D = new double[1][];
        if (numPadded.length > 0) {
            D[0] = new double[]{numPadded[0]};
        } else {
            /*
                We don't assign it an empty array because this system
                is not 'null'. It just doesn't have a non-zero D
                matrix. Thus, it should have a non-zero shape so that
                it can be operated on by functions like 'ss2tf'
             */
            D[0] = new double[]{0};
        }
        int k = den.length;
        if (k == 1) {
            return new StateSpace(new double[][]{{0.0}}, new double[][]{{0.0}}, new double[][]{{0.0}}, D);
        }
        double[] fRow = new double[k - 1];
        System.arraycopy(den, 1, fRow, 0, fRow.length);
        NumArrays.multiplyElementWiseInPlace(fRow, -1.0);

        double[][] eye = Matrices.identity(k - 2, k - 1).getAs2DArray();
        double[][] A = new double[eye.length + 1][];
        A[0] = fRow;
        for (int i = 0; i < eye.length; ++i) {
            A[i + 1] = Arrays.copyOf(eye[i], eye[0].length);
        }
        double[][] B = Matrices.identity(k - 1, 1).getAs2DArray();
        double[][] C = new double[1][];
        double[][] outer = NumArrays.outer(new double[]{numPadded[0]}, Arrays.copyOfRange(den, 1, den.length));
        C[0] = NumArrays.subtract(Arrays.copyOfRange(numPadded, 1, numPadded.length), outer[0]);
        return new StateSpace(A, B, C, D);
    }

    public void normalize() {
        this.rf.normalize();
    }
    // TODO add getOrder

    @Override
    public TransferFunction toTransferFunction() {
        return this;
    }


    public SISOTimeResponse simulateTimeResponse(double[] input, double[] time) {
        return simulateTimeResponse(input, time, IntegrationMethod.INTERPOLATION);
    }

    public SISOTimeResponse simulateTimeResponse(double[] input, double[] time,
                                                 double[] initialConditions) {
        return simulateTimeResponse(input, time, initialConditions, IntegrationMethod.INTERPOLATION);
    }

    public SISOTimeResponse simulateTimeResponse(double[] input, double[] time,
                                                 IntegrationMethod integrationMethod) {
        return simulateTimeResponse(input, time, null, integrationMethod);
    }

    public SISOTimeResponse simulateTimeResponse(double[] input, double[] time,
                                                 double[] initialConditions,
                                                 IntegrationMethod integrationMethod) {
        double[][] U = new double[1][time.length];
        U[0] = input;
        TimeResponse tr = lsim(U, time, initialConditions, this.toStateSpace(), integrationMethod);
        return new SISOTimeResponse(tr.getTime(), tr.getResponse()[0], tr.getEvolutionOfStateVector());
    }
}
