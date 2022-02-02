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

    /**
     * The {@code Frequency Response} holds the {@link Complex }results of the evaluation of the
     * {@link TransferFunction} at a given set of frequencies.
     */
    public static class FrequencyResponse {
        private Complex[] response;
        private double[] w;

        FrequencyResponse(Complex[] response, double[] w) {
            this.response = response;
            this.w = w;
        }

        /**
         * Complex response of the system.
         * @return The {@link Complex} response of the system.
         */
        public Complex[] getResponse() {
            return response;
        }

        /**
         * Frequencies at which the system was evaluated.
         * @return An array of frequencies at which the system was evaluated.
         */
        public double[] getFrequencies() {
            return w;
        }
    }

    /**
     * The {@code BodeResponse} class holds the magnitude and frequency response after evaluating the given
     * {@link TransferFunction} at an array of frequencies.
     */
    public static class BodeResponse {
        private double[] magnitudeIndB;
        private double[] phase;
        private double[] w;

        BodeResponse(double[] magnitudeIndB, double[] phase, double[] w) {
            this.magnitudeIndB = magnitudeIndB;
            this.phase = phase;
            this.w = w;
        }

        /**
         * Magnitude response of the system.
         * @return The magnitude response of the system. <br>
         * {@code Mag<sub>dB</sub> = 20 * log10(abs(response))}.
         */
        public double[] getMagnitudeIndB() {
            return magnitudeIndB;
        }

        /**
         * Phase response of the system.
         * @return The wrapped phase response of the system is degrees. The wrapped phase only goes from -180° to 180°. <br>
         * To unwrap the phase the {@link TransferFunction#unwrapPhase(double[])} can be used.
         */
        public double[] getPhaseInDegrees() {
            return phase;
        }

        /**
         * Frequencies at which the system was evaluated.
         * @return An array of frequencies at which the system was evaluated.
         */
        public double[] getFrequencies() {
            return w;
        }
    }

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

    public Complex[] getZeros() {
        return rf.getZeros();
    }

    public Complex[] getPoles() {
        return rf.getPoles();
    }

    public Polynomial getNumerator() {
        return rf.getNumerator();
    }

    public Polynomial getDenominator() {
        return rf.getDenominator();
    }

    public Complex evaluateAt(double w) {
        return rf.evaluateAt(0.0, w);
    }

    public double[] getMagnitudeAt(double[] w) {
        double[] magnitude = new double[w.length];
        for (int i = 0; i < w.length; ++i) {
            magnitude[i] = this.evaluateAt(w[i]).abs();
        }
        return magnitude;
    }

    /***
     * Calculate the system wrapped phase response.
     * @param frequencies the frequencies where the phase needs to be calculated at.
     * @return The phase response of the system in rad / s.
     */
    public double[] getPhaseAt(double[] frequencies) {
        double[] phase = new double[frequencies.length];
        for (int i = 0; i < frequencies.length; ++i) {
            phase[i] = this.evaluateAt(frequencies[i]).arg();
        }
        return phase;
    }

    /***
     * Calculate the system wrapped phase response.
     * @param frequencies the frequencies where the phase needs to be calculated at.
     * @return The phase response of the system in degrees.
     */
    public double[] getPhaseInDegreesAt(double[] frequencies) {
        double[] phase = new double[frequencies.length];
        for (int i = 0; i < frequencies.length; ++i) {
            phase[i] = this.evaluateAt(frequencies[i]).arg() * (180 / Math.PI);
        }
        return phase;
    }

    public FrequencyResponse getFrequencyResponse() {
        double[] frequencies = findFrequencies(200);
        return this.getFrequencyResponse(frequencies);
    }

    public FrequencyResponse getFrequencyResponse(int numberOfPoints) {
        double[] frequencies = findFrequencies(numberOfPoints);
        return this.getFrequencyResponse(frequencies);
    }

    public FrequencyResponse getFrequencyResponse(double[] w) {
        Complex[] response = new Complex[w.length];
        for (int i = 0; i < w.length; ++i) {
            response[i] = this.evaluateAt(w[i]);
        }
        return new FrequencyResponse(response, w);
    }

    public BodeResponse getBode() {
        double[] frequencies = findFrequencies(200);
        return this.getBode(frequencies);
    }

    public BodeResponse getBode(int numberOfPoints) {
        double[] frequencies = findFrequencies(numberOfPoints);
        return this.getBode(frequencies);
    }

    public BodeResponse getBode(double[] w) {
        double[] magnitudeIndB = new double[w.length];
        double[] phaseInDegrees = new double[w.length];
        for (int i = 0; i < w.length; ++i) {
            magnitudeIndB[i] = 20 * Math.log10(this.getMagnitudeAt(w[i]));
            phaseInDegrees[i] = this.evaluateAt(w[i]).arg() * (180 / Math.PI);
        }
        return new BodeResponse(magnitudeIndB, phaseInDegrees, w);
    }

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

    public TransferFunction add(final TransferFunction tf) {
        return new TransferFunction(rf.add(tf.rf));
    }

    public TransferFunction add(double d) {
        return new TransferFunction(rf.add(d));
    }

    public TransferFunction subtract(TransferFunction tf) {
        return new TransferFunction(rf.subtract(tf.rf));
    }

    public TransferFunction multiply(TransferFunction tf) {
        return new TransferFunction(rf.multiply(tf.rf));
    }

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

    private static Polynomial getPolynomialMagnitude(Polynomial polynomial) {
        double[] num = polynomial.getCoefficients();
        Complex[] numjw = evaluateAtjw(num);
        return new Polynomial(norm(numjw));
    }

    /***
     * Calculates the squared magnitude element wise
     * @param a
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

    private static int getComplexSign(int power) {
        int result = power - 4 * (power / 4);
        return result == 2 || result == 3 ? -1 : 1;
    }

    public double[] getAllGainCrossoverFrequencies() {
        return this.getAllGainCrossoverFrequencies(Math.sqrt(ConstantsETK.DOUBLE_EPS));
    }


    // 		  |Num(s)|
    // H(s) = -------- = 1
    //        |Den(s)|
    //
    // |Num(s)| - |Den(s)| = 0
    //
    // Find roots of the resultant polynomial

    public double[] getAllGainCrossoverFrequencies(double tol) {

        Polynomial magPolyNum = getPolynomialMagnitude(rf.getNumerator());
        Polynomial magPolyDen = getPolynomialMagnitude(rf.getDenominator());

        Complex[] solution = magPolyNum.subtract(magPolyDen).getRoots();
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

    public double getGainCrossoverFrequency() {
        double[] freqs = this.getAllGainCrossoverFrequencies();
        // Get the first one were the magnitude was decreasing
        for (double f : freqs) {
            double slope = this.getMagnitudeAt(f) - this.getMagnitudeAt(f + 1e-12);
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
    public double getPhaseCrossoverFrequency() {
        return this.getPhaseCrossoverFrequency(Math.sqrt(ConstantsETK.DOUBLE_EPS));
    }

    /**
     * Phase crossover frequency.
     * @return The first frequency at which the phase crosses -180°.
     */
    public double getPhaseCrossoverFrequency(double tol) {
        double[] frequencies = this.getAllPhaseCrossoverFrequencies(tol);
        // Get the first one where the magnitude was decreasing
        for (double freq : frequencies) {
            double slope = this.getMagnitudeAt(freq) - this.getMagnitudeAt(freq + 1e-12);
            if (slope > 0.0) {
                return freq;
            }
        }
        return Double.NaN;
    }

    public double[] getAllPhaseCrossoverFrequencies() {
        return this.getAllPhaseCrossoverFrequencies(Math.sqrt(ConstantsETK.DOUBLE_EPS));
    }

    /**
     * Finds all the phase crossover frequencies by solving H(s) = H(-s). The phase crossover frequencies are defined
     * as the frequencies at which the phase of the system crosses -180°. Typically, the first frequency at which
     * the phase crosses -180° is the frequency of interest and in that case
     * {@link TransferFunction#getPhaseCrossoverFrequency()} can be used. <br>
     *
     * <pre>
     * Let's define our system as:
     *
     *         Num(jw)
     * H(jw) = -------
     *         Den(jw)
     *
     * The goal is to solve
     *
     *         Num(jw) * Den'(jw)
     * H(jw) = ------------------
     *         Den(jw) * Den'(jw)
     *
     * Where H'(jw) = Conjugate of H(jw).
     *
     * After multiplying the denominator with its conjugate, it becomes a real
     * number and it won't affect the phase, therefore it doesn't have to be taken
     * into account from this point on. We can also drop Num'(jw) * Den'(jw) since
     * it's there to zero out the real part of Num(jw) * Den'(jw) - Num'(jw) * Den'(jw),
     * which can be also accomplished by setting the real part of Num(jw) * Den'(jw)
     * to zero.
     *
     * Let Nump(jw) = Num(jw) * Den'(jw) thus
     *     Nump(jw) = Real(Nump(jw)) + j * Imag(Nump(jw))
     *
     * We can equate the imaginary part of Num'(jw) to 0 and solve for w
     *
     * The roots of Imag(Nump(jw)) are the solution vector containing the gain
     * crossover frequencies.
     * </pre>
     * @return An array of the phase crossover frequencies in ascending order.
     */
    public double[] getAllPhaseCrossoverFrequencies(double tol) {
        Complex[] num = evaluateAtjw(rf.getNumerator().getCoefficients());
        Complex[] den = evaluateAtjw(rf.getDenominator().getCoefficients());

        double[] convLHS = NumArrays.convolution(ComplexArrays.imag(num), ComplexArrays.real(den));
        double[] convRHS = NumArrays.convolution(ComplexArrays.real(num), ComplexArrays.imag(den));

        Polynomial a = new Polynomial(convLHS);
        Polynomial b = new Polynomial(convRHS);
        a.subtractEquals(b);

        Complex[] roots = a.getRoots();
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

    public double getMagnitudeAt(double f) {
        return this.evaluateAt(f).abs();
    }

    public Margins getMargins() {
        double wcg = this.getGainCrossoverFrequency();
        double wcp = this.getPhaseCrossoverFrequency();

        double pm = Double.POSITIVE_INFINITY;
        if (!Double.isNaN(wcg)) {
            double phase = this.getPhaseInDegreesAt(wcg);
            pm = 180.0 + phase;
        }
        double gm = Double.POSITIVE_INFINITY;
        if (!Double.isNaN(wcp)) {
            double mag = this.getMagnitudeAt(wcp);
            gm = 1.0 / mag;
        }

        Margins m = new Margins(gm, pm, wcp, wcg);
        return m;
    }

    public void substituteInPlace(double d) {
        rf.substituteInPlace(d);
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
     * @param f
     * @return the phase of the system in degrees
     */
    public double getPhaseInDegreesAt(double f) {
        Complex[] zeros = this.getZeros();
        Complex[] poles = this.getPoles();

        double phase = calculatePhase(zeros, f);
        phase -= calculatePhase(poles, f);

        return Math.toDegrees(phase);
    }


    /***
     *
     * @param roots
     * @param f
     * @return
     */
    private static double calculatePhase(Complex[] roots, double f) {
        double phase = 0.0;
        for (int i = 0; i < roots.length; ++i) {
            if (roots[i].imag() == 0) {
                // Single real root
                double alpha = -roots[i].real();
                phase += Math.atan2(f, alpha);
            } else {
                // Complex pair of roots
                double alpha = -roots[i].real();
                double beta = roots[i].imag();
                double fn = MathETK.hypot(alpha, beta);
                double zeta = alpha / fn;
                double ratio = f / fn;
                phase += Math.atan2(2.0 * zeta * ratio, 1 - ratio * ratio);
                // Roots are in pairs. Skip the conjugate since it was
                // already taken into account in the calculation above
                ++i;
            }
        }
        return phase;
    }

    public boolean isProper() {
        return rf.isProper();
    }

    public boolean isStrictlyProper() {
        return rf.isStrictlyProper();
    }

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
            throw new ImproperTransferFunctionException("Transfer function must be proper");
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

        double[][] eye = Matrices.Identity(k - 2, k - 1).getAs2DArray();
        double[][] A = new double[eye.length + 1][];
        A[0] = fRow;
        for (int i = 0; i < eye.length; ++i) {
            A[i + 1] = Arrays.copyOf(eye[i], eye[0].length);
        }
        double[][] B = Matrices.Identity(k - 1, 1).getAs2DArray();
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


    public SingleInputSingleOutputTimeResponse simulateTimeResponse(double[] input, double[] time) {
        return simulateTimeResponse(input, time, IntegrationMethod.INTERPOLATION);
    }

    public SingleInputSingleOutputTimeResponse simulateTimeResponse(double[] input, double[] time,
                                                                    double[] initialConditions) {
        return simulateTimeResponse(input, time, initialConditions, IntegrationMethod.INTERPOLATION);
    }

    public SingleInputSingleOutputTimeResponse simulateTimeResponse(double[] input, double[] time,
                                                                    IntegrationMethod integrationMethod) {
        return simulateTimeResponse(input, time, null, integrationMethod);
    }

    public SingleInputSingleOutputTimeResponse simulateTimeResponse(double[] input, double[] time,
                                                                    double[] initialConditions,
                                                                    IntegrationMethod integrationMethod) {
        double[][] U = new double[1][time.length];
        U[0] = input;
        TimeResponse tr = lsim(U, time, initialConditions, this.toStateSpace(), integrationMethod);
        return new SingleInputSingleOutputTimeResponse(tr.getTime(), tr.getResponse()[0], tr.getEvolutionOfStateVector());
    }

    public static void main(String[] args) {
//        System.out.println(roundFrequency(1.17));
//        System.out.println(roundFrequency(-0.23));
//        System.out.println(roundFrequency(-0.5));
//        System.out.println(roundFrequency(-0.57));
//        System.out.println(roundFrequency(-0.51));
//        System.out.println(roundFrequency(2.5));
//        System.out.println(roundFrequency(2.51));
//        System.out.println(roundFrequency(2.49));
//        System.out.println(roundFrequency(3.5));
//        System.out.println(roundFrequency(3.49));
//        TransferFunction tf1 = new TransferFunction(new double[]{1000, 0}, new double[]{1, 25, 100, 9, 4});
//        System.out.println(tf1);
//        System.out.println(tf1.getMargins());
//
//        double phase = tf1.getPhaseInDegreesAt(70.4);
//        System.out.println(phase);
//
//        TransferFunction tf2 = new TransferFunction(new double[]{1, 3, 3}, new double[]{1, 2, 1});
//        System.out.println(tf2.toStateSpace());
//
//        TransferFunction tf3 = new TransferFunction(new double[]{5, 3, 4}, new double[]{8, 2, 9, 10});
//        System.out.println(tf3.toStateSpace());
//
//        TransferFunction tf4 = new TransferFunction(new double[]{5}, new double[]{3});
//        System.out.println(tf4.toStateSpace());
//
//        TransferFunction tf5 = new TransferFunction(new double[]{1.0, 3, 3}, new double[]{1.0, 2.0, 1});
//        tf5.step();
//
//        double[] timePoints = {0.0, 0.0707070707070707, 0.1414141414141414, 0.2121212121212121, 0.2828282828282828,
//                0.35353535353535354, 0.4242424242424242, 0.4949494949494949, 0.5656565656565656, 0.6363636363636364,
//                0.7070707070707071, 0.7777777777777778, 0.8484848484848484, 0.9191919191919191, 0.9898989898989898,
//                1.0606060606060606, 1.1313131313131313, 1.202020202020202, 1.2727272727272727, 1.3434343434343434,
//                1.4141414141414141, 1.4848484848484849, 1.5555555555555556, 1.6262626262626263, 1.6969696969696968,
//                1.7676767676767675, 1.8383838383838382, 1.909090909090909, 1.9797979797979797, 2.0505050505050506,
//                2.121212121212121, 2.191919191919192, 2.2626262626262625, 2.333333333333333, 2.404040404040404,
//                2.4747474747474745, 2.5454545454545454, 2.616161616161616, 2.686868686868687, 2.7575757575757573,
//                2.8282828282828283, 2.898989898989899, 2.9696969696969697, 3.04040404040404, 3.111111111111111,
//                3.1818181818181817, 3.2525252525252526, 3.323232323232323, 3.3939393939393936, 3.4646464646464645,
//                3.535353535353535, 3.606060606060606, 3.6767676767676765, 3.7474747474747474, 3.818181818181818,
//                3.888888888888889, 3.9595959595959593, 4.03030303030303, 4.101010101010101, 4.171717171717171,
//                4.242424242424242, 4.313131313131313, 4.383838383838384, 4.454545454545454, 4.525252525252525,
//                4.595959595959596, 4.666666666666666, 4.737373737373737, 4.808080808080808, 4.878787878787879,
//                4.949494949494949, 5.02020202020202, 5.090909090909091, 5.161616161616162, 5.232323232323232,
//                5.303030303030303, 5.373737373737374, 5.444444444444445, 5.515151515151515, 5.585858585858586,
//                5.656565656565657, 5.727272727272727, 5.797979797979798, 5.8686868686868685, 5.9393939393939394,
//                6.0101010101010095, 6.08080808080808, 6.151515151515151, 6.222222222222222, 6.292929292929292,
//                6.363636363636363, 6.434343434343434, 6.505050505050505, 6.575757575757575, 6.646464646464646,
//                6.717171717171717, 6.787878787878787, 6.858585858585858, 6.929292929292929, 7.0};
//
//        TransferFunction tf6 = new TransferFunction(new double[]{1.0, 3.0, 3.0}, new double[]{1.0, 2.0, 1.0});
//        double[] yOut = tf6.step(timePoints, new double[]{1.0, 0.0}).getResponse();
//
//        System.out.println(Arrays.toString(yOut));
//
//        TransferFunction tf7 = new TransferFunction(new double[]{1.0, 3.0, 3.0}, new double[]{1.0, 2.0, 1.0});
//        double[] yOut2 = tf7.step(new double[]{0.0}).getResponse();
//
//        System.out.println(Arrays.toString(yOut2));
//
//        TransferFunction tf8 = new TransferFunction(new double[] {1, 1, 0}, new double[] {1, 8, 25, 12, 1, 9 , 8, 10});
//        double[] wn = tf8.findFrequencies(9);
//
//        System.out.println(Arrays.toString(wn));
//
//        TransferFunction tf9 = new TransferFunction(new double[] {1, 0.1, 7.5}, new double[] {1, 0.12, 9, 0, 0});
//
//        System.out.println(Arrays.toString(tf9.findFrequencies(9)));
//
//        // Specs for band pass filter
//        FilterSpecs.BandPassSpecs bpSpecs = new FilterSpecs.BandPassSpecs();
//        // The bandwidth of the filter starts at the LowerPassBandFrequency and
//        // ends at the UpperPassBandFrequency. The filter has lower stop band
//        // which is set LowerStopBandFrequency and the upper stop band can be set
//        // with UpperStopBandFrequency. The attenuation at the stop bands can be
//        // set with the LowerStopBandAttenuation and UpperStopBandAttenuation
//        // respectively. In a frequency spectrum, the order of the frequencies will be:
//        // LowerStopBandFrequency < LowerPassBandFrequency < UpperPassBandFrequency <
//        // UpperStopBandFrequency
//        bpSpecs.setLowerPassBandFrequency(190.0); // 190 Hz lower pass band frequency
//        bpSpecs.setUpperPassBandFrequency(210.0); // 210 Hz upper pass band frequency
//        bpSpecs.setLowerStopBandFrequency(180.0); // 180 Hz lower stop band frequency
//        bpSpecs.setUpperStopBandFrequency(220.0); // 220 Hz upper stop band frequency
//        bpSpecs.setPassBandRipple(0.2); // 0.2 dB gain/ripple refer to note
//        bpSpecs.setStopBandAttenuation(20.0); // 20 dB attenuation in the stop band
//
//        FilterOrderResults.OrderAndCutoffFrequencies nW0W1 = Elliptic.ellipord(bpSpecs);
//        TransferFunction el = Elliptic.newBandPass(nW0W1.getOrder(), bpSpecs.getPassBandRipple(),
//                bpSpecs.getStopBandAttenuation(), nW0W1.getLowerCutoffFrequency(), nW0W1.getUpperCutoffFrequency());
//
//        System.out.println(Arrays.toString(el.findFrequencies(9)));
//
//        FilterSpecs.BandStopSpecs bsSpecs = new FilterSpecs.BandStopSpecs();
//        // The notch of the filter starts at the LowerStopBandFrequency and
//        // ends at the UpperStopBandFrequency. The filter has lower pass band
//        // which is set LowerPassBandFrequency and the upper pass band can be set
//        // with UpperPassBandFrequency. The attenuation at the notch can be
//        // set with the StopBandAttenuation parameter and the attenuation/ripple
//        // in the pass band can be set with the PassBandRipple parameter.
//        // In a frequency spectrum, the order of the frequencies will be:
//        // LowerPassBandFrequency < LowerStopBandFrequency < UpperStopBandFrequency <
//        // UpperPassBandFrequency
//        bsSpecs.setLowerPassBandFrequency(3.6e3); // 3600 Hz lower pass band frequency
//        bsSpecs.setUpperPassBandFrequency(9.1e3); // 9100 Hz lower pass band frequency
//        bsSpecs.setLowerStopBandFrequency(5.45e3); // 5450 Hz lower stop band frequency
//        bsSpecs.setUpperStopBandFrequency(5.90e3); // 5900 Hz upper stop band frequency
//        bsSpecs.setPassBandRipple(0.5); // 1.5 dB gain/ripple refer to note
//        bsSpecs.setStopBandAttenuation(38.0); // 38 db attenuation at the notch
//
//        nW0W1 = Elliptic.ellipord(bsSpecs);
//        el = Elliptic.newBandStop(nW0W1.getOrder(), bsSpecs.getPassBandRipple(),
//                bsSpecs.getStopBandAttenuation(), nW0W1.getLowerCutoffFrequency(), nW0W1.getUpperCutoffFrequency());
//
//        System.out.println(Arrays.toString(el.findFrequencies(9)));
//
//        TransferFunction tf10 = new TransferFunction(new double[] {1, 0.1, 7.5, 0, 1, 0, 0 ,-10},
//                new double[] {1, 0, -20, 1e13, 0, 8, 1, 0.12, 9, 0, 0});
//
//        System.out.println(Arrays.toString(tf10.findFrequencies(9)));

        System.out.println(getComplexSign(2));
//		double[] phase = tf1.getPhaseAt(logspace);
//		System.out.println(Arrays.toString(phase));
//		for(int i = 0; i < phase.length; ++i) {
//			System.out.printf("%.4f %.4f%n", logspace[i], phase[i]);
//		}
//		
//		System.out.println();
//		for (double val : tf1.getAllPhaseCrossoverFrequencies()) {
//			System.out.print(val + " ");
//		}
//		System.out.println();
//		TransferFunction tf2 = new TransferFunction(new double[] { 12, 6 }, new double[] { 1, 5, 10 });
//		for (double val : tf2.getAllGainCrossoverFrequencies()) {
//			System.out.print(val + " ");
//		}
//		System.out.println();
//		for (double val : tf2.getAllPhaseCrossoverFrequencies()) {
//			System.out.print(val + " ");
//		}
//		System.out.println();
//		System.out.println(tf1.getMargins());
//		TransferFunction tf = new TransferFunction(zeros.toArray(new Complex[zeros.size()]),
//				poles.toArray(new Complex[poles.size()]));
//
//		double x = tf.getMagnitudeAt(1.0);
//		double y = tf.getPhaseAt(1.0);
//		System.out.println(x);
//		System.out.println(y);
//
//		System.out.println(tf1.toString());
//
//		double[] freq = ArrayUtils.logspace(-3, 3, 100);
//
//		freq = ArrayUtils.logspace(-3, 3, 10000000);

    }
}
