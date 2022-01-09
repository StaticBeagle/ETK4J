package com.wildbitsfoundry.etk4j.control;

import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrices;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;
import com.wildbitsfoundry.etk4j.signals.filters.Elliptic;
import com.wildbitsfoundry.etk4j.signals.filters.FilterOrderResults;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;
import com.wildbitsfoundry.etk4j.util.NumArrays;

import java.util.Arrays;

/***
 *
 * @author StaticBeagle
 *
 */
public class TransferFunction extends LinearTimeInvariantSystem {

    public static class FrequencyResponse {
        private Complex[] response;
        private double[] frequencies;

        // TODO rename all frequencies to wn?
        FrequencyResponse(Complex[] response, double[] frequencies) {
            this.response = response;
            this.frequencies = frequencies;
        }

        /***
         * Returns the complex magnitude of a system
         * @return
         */
        public Complex[] getResponse() {
            return response;
        }
        // TODO document this
        public double[] getFrequencies() {
            return frequencies;
        }
    }

    public static class BodeResponse {
        private double[] magnitudeIndB;
        private double[] phase;
        private double[] frequencies;

        // TODO rename all frequencies to wn?
        BodeResponse(double[] magnitudeIndB, double[] phase, double[] frequencies) {
            this.magnitudeIndB = magnitudeIndB;
            this.phase = phase;
            this.frequencies = frequencies;
        }

        /***
         * Returns the complex magnitude of a system
         * @return
         */
        public double[] getMagnitudeIndB() {
            return magnitudeIndB;
        }

        public double[] getPhaseInDegrees() {
            return phase;
        }
        // TODO document this
        public double[] getFrequencies() {
            return frequencies;
        }
    }

    private RationalFunction _rf;

    /***
     * Constructs a transfer function from poles and zeros
     * @param zeros - zeros of the transfer function
     * @param poles - poles of the transfer function
     */
    public TransferFunction(Complex[] zeros, Complex[] poles) {
        _rf = new RationalFunction(zeros, poles);
    }

    public TransferFunction(double gain, Complex[] poles) {
        Polynomial num = new Polynomial(gain);
        Polynomial den = new Polynomial(poles);

        _rf = new RationalFunction(num, den);
    }

    public TransferFunction(TransferFunction tf) {
        _rf = new RationalFunction(tf._rf);
    }

    public TransferFunction(Polynomial numerator, Polynomial denominator) {
        _rf = new RationalFunction(numerator, denominator);
    }

    public TransferFunction(double[] numerator, double[] denominator) {
        _rf = new RationalFunction(numerator, denominator);
    }

    public TransferFunction(RationalFunction rf) {
        _rf = rf;
    }

    public TransferFunction(ZeroPoleGain zpk) {
        this(zpk.getZeros(), zpk.getPoles(), zpk.getGain());
    }

    public TransferFunction(Complex[] zeros, Complex[] poles, double gain) {
        _rf = new RationalFunction(zeros, poles, gain);
    }

    @Override
    public ZeroPoleGain toZeroPoleGain() {
        Complex[] zeros = _rf.getZeros();
        Complex[] poles = _rf.getPoles();
        double k = _rf.getNumerator().getCoefficientAt(0);
        return new ZeroPoleGain(zeros, poles, k);
    }

    public Complex[] getZeros() {
        return _rf.getZeros();
    }

    public Complex[] getPoles() {
        return _rf.getPoles();
    }

    public Polynomial getNumerator() {
        return _rf.getNumerator();
    }

    public Polynomial getDenominator() {
        return _rf.getDenominator();
    }

    public Complex evaluateAt(double f) {
        return _rf.evaluateAt(0.0, f);
    }

    public double[] getMagnitudeAt(double[] frequencies) {
        double[] magnitude = new double[frequencies.length];
        for (int i = 0; i < frequencies.length; ++i) {
            magnitude[i] = this.evaluateAt(frequencies[i]).abs();
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

    // TODO implement getFrequencyResponse(void)?
    public FrequencyResponse getFrequencyResponse() {
        double[] frequencies = findFrequencies(200);
        return this.getFrequencyResponse(frequencies);
    }

    public FrequencyResponse getFrequencyResponse(int numberOfPoints) {
        double[] frequencies = findFrequencies(numberOfPoints);
        return this.getFrequencyResponse(frequencies);
    }

    // TODO rename frequencies to wn?
    public FrequencyResponse getFrequencyResponse(double[] frequencies) {
        Complex[] response = new Complex[frequencies.length];
        for (int i = 0; i < frequencies.length; ++i) {
            response[i] = this.evaluateAt(frequencies[i]);
        }
        return new FrequencyResponse(response, frequencies);
    }

    public BodeResponse getBode() {
        double[] frequencies = findFrequencies(200);
        return this.getBode(frequencies);
    }

    public BodeResponse getBode(int numberOfPoints) {
        double[] frequencies = findFrequencies(numberOfPoints);
        return this.getBode(frequencies);
    }

    // TODO rename frequencies to wn?
    public BodeResponse getBode(double[] frequencies) {
        double[] magnitudeIndB = new double[frequencies.length];
        double[] phaseInDegrees = new double[frequencies.length];
        for (int i = 0; i < frequencies.length; ++i) {
            magnitudeIndB[i] = 20 * Math.log10(this.getMagnitudeAt(frequencies[i]));
            phaseInDegrees[i] = this.evaluateAt(frequencies[i]).arg() * (180 / Math.PI);
        }
        return new BodeResponse(magnitudeIndB, phaseInDegrees, frequencies);
    }

    // TODO rename to numberOfPoints?
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
        return new TransferFunction(_rf.add(tf._rf));
    }

    public TransferFunction add(double d) {
        return new TransferFunction(_rf.add(d));
    }

    public TransferFunction subtract(TransferFunction tf) {
        return new TransferFunction(_rf.subtract(tf._rf));
    }

    public TransferFunction multiply(TransferFunction tf) {
        return new TransferFunction(_rf.multiply(tf._rf));
    }

    public TransferFunction multiply(double d) {
        return new TransferFunction(_rf.multiply(d));
    }

    @Override
    public String toString() {
        String num = _rf.getNumerator().toString();
        String den = _rf.getDenominator().toString();

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
        Complex[] numjw = evalAtjw(num);
        return new Polynomial(norm(numjw));
    }

    /***
     * Calculates the squared magnitude element wise
     * @param a
     * @return [norm0, norm1, .... normn]
     */
    private static double[] norm(Complex[] a) {
        Complex[] mag = ComplexArrays.convolution(a, ComplexArrays.conj(a));
        double[] coefs = new double[mag.length];
        for (int i = 0; i < mag.length; ++i) {
            coefs[i] = mag[i].real();
        }
        return coefs;
    }

    private static Complex[] evalAtjw(double[] poly) {
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


    // 		  |Num(s)|
    // H(s) = -------- = 1
    //        |Den(s)|
    //
    // |Num(s)| - |Den(s)| = 0
    //
    // Find roots of the resultant polynomial

    public double[] getAllGainCrossoverFrequencies() {

        Polynomial magPolyNum = getPolynomialMagnitude(_rf.getNumerator());
        Polynomial magPolyDen = getPolynomialMagnitude(_rf.getDenominator());

        Complex[] solution = magPolyNum.subtract(magPolyDen).getRoots();
        double[] wgc = new double[solution.length];
        int j = 0;
        for (int i = 0; i < solution.length; i++) {
            double real = solution[i].real();
            if (solution[i].imag() == 0 && real > 0) {
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
        return Double.POSITIVE_INFINITY;
    }


    /***
     * Finds all the phase cross over frequencies by solving H(s) = H(-s)
     *
     * <pre>
     * This is equivalent to:
     *
     *         Num(jw) * Den'(jw) - Num'(jw) * Den'(jw)
     * H(jw) = ----------------------------------------
     *                    Den(jw) * Den'(jw)
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
     * @return Returns an array of values that are the phase crossover
     *         frequencies
     */
    public double[] getAllPhaseCrossoverFrequencies() {
        Complex[] num = evalAtjw(_rf.getNumerator().getCoefficients());
        Complex[] conjDen = ComplexArrays.conj(evalAtjw(_rf.getDenominator().getCoefficients()));

        Complex[] conv = ComplexArrays.convolution(num, conjDen);
        double[] imag = new double[conv.length];
        for (int i = 0; i < conv.length; ++i) {
            imag[i] = conv[i].imag();
        }

        Complex[] solution = new Polynomial(imag).getRoots();

        double[] wpc = new double[solution.length];
        int j = 0;
        for (int i = 0; i < solution.length; i++) {
            double real = solution[i].real();
            if (solution[i].imag() == 0 && real > 0) {
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

    public double getPhaseCrossoverFrequency() {
        double[] freqs = this.getAllPhaseCrossoverFrequencies();
        // Get the first one were the magnitude was decreasing
        for (double f : freqs) {
            double slope = this.getMagnitudeAt(f) - this.getMagnitudeAt(f + 1e-12);
            if (slope > 0.0) {
                return f;
            }
        }
        return Double.POSITIVE_INFINITY;
    }

    public Margins getMargins() {
        double wcg = this.getGainCrossoverFrequency();
        double wcp = this.getPhaseCrossoverFrequency();

        double pm = wcg;
        if (wcg != Double.POSITIVE_INFINITY) {
            double phase = this.getPhaseInDegreesAt(wcg);
            pm = 180.0 + phase;
        }
        double gm = wcp;
        if (wcp != Double.POSITIVE_INFINITY) {
            double mag = this.getMagnitudeAt(wcp);
            gm = 20 * Math.log10(1.0 / mag);
        }

        Margins m = new Margins(gm, pm, wcg, wcp);
        return m;
    }

    public void substituteInPlace(double d) {
        _rf.substituteInPlace(d);
    }

    /***
     * Calculates the phase at of the system a given frequency. </br>
     * This operation uses the zeros and poles of the system to calculate the phase as: </br>
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
        return _rf.isProper();
    }

    public boolean isStrictlyProper() {
        return _rf.isStrictlyProper();
    }

    /***
     * Transform SISO only single input single output TFs
     * @return
     */
    @Override
    public StateSpace toStateSpace() {

        TransferFunction tf = new TransferFunction(this);

        tf.normalize();
        if (!tf.isProper()) {
            throw new ImproperTransferFunctionException();
        }
        double[] num = tf._rf.getNumerator().getCoefficients();
        double[] den = tf._rf.getDenominator().getCoefficients();

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
        NumArrays.multiplyInPlace(fRow, -1.0);

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
        this._rf.normalize();
    }
    // TODO add getOrder

    // TODO this should be public to comply with Liskov
    @Override
    protected TransferFunction toTransferFunction() {
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
        System.out.println(roundFrequency(1.17));
        System.out.println(roundFrequency(-0.23));
        System.out.println(roundFrequency(-0.5));
        System.out.println(roundFrequency(-0.57));
        System.out.println(roundFrequency(-0.51));
        System.out.println(roundFrequency(2.5));
        System.out.println(roundFrequency(2.51));
        System.out.println(roundFrequency(2.49));
        System.out.println(roundFrequency(3.5));
        System.out.println(roundFrequency(3.49));
        TransferFunction tf1 = new TransferFunction(new double[]{1000, 0}, new double[]{1, 25, 100, 9, 4});
        System.out.println(tf1);
        System.out.println(tf1.getMargins());

        double phase = tf1.getPhaseInDegreesAt(70.4);
        System.out.println(phase);

        TransferFunction tf2 = new TransferFunction(new double[]{1, 3, 3}, new double[]{1, 2, 1});
        System.out.println(tf2.toStateSpace());

        TransferFunction tf3 = new TransferFunction(new double[]{5, 3, 4}, new double[]{8, 2, 9, 10});
        System.out.println(tf3.toStateSpace());

        TransferFunction tf4 = new TransferFunction(new double[]{5}, new double[]{3});
        System.out.println(tf4.toStateSpace());

        TransferFunction tf5 = new TransferFunction(new double[]{1.0, 3, 3}, new double[]{1.0, 2.0, 1});
        tf5.step();

        double[] timePoints = {0.0, 0.0707070707070707, 0.1414141414141414, 0.2121212121212121, 0.2828282828282828,
                0.35353535353535354, 0.4242424242424242, 0.4949494949494949, 0.5656565656565656, 0.6363636363636364,
                0.7070707070707071, 0.7777777777777778, 0.8484848484848484, 0.9191919191919191, 0.9898989898989898,
                1.0606060606060606, 1.1313131313131313, 1.202020202020202, 1.2727272727272727, 1.3434343434343434,
                1.4141414141414141, 1.4848484848484849, 1.5555555555555556, 1.6262626262626263, 1.6969696969696968,
                1.7676767676767675, 1.8383838383838382, 1.909090909090909, 1.9797979797979797, 2.0505050505050506,
                2.121212121212121, 2.191919191919192, 2.2626262626262625, 2.333333333333333, 2.404040404040404,
                2.4747474747474745, 2.5454545454545454, 2.616161616161616, 2.686868686868687, 2.7575757575757573,
                2.8282828282828283, 2.898989898989899, 2.9696969696969697, 3.04040404040404, 3.111111111111111,
                3.1818181818181817, 3.2525252525252526, 3.323232323232323, 3.3939393939393936, 3.4646464646464645,
                3.535353535353535, 3.606060606060606, 3.6767676767676765, 3.7474747474747474, 3.818181818181818,
                3.888888888888889, 3.9595959595959593, 4.03030303030303, 4.101010101010101, 4.171717171717171,
                4.242424242424242, 4.313131313131313, 4.383838383838384, 4.454545454545454, 4.525252525252525,
                4.595959595959596, 4.666666666666666, 4.737373737373737, 4.808080808080808, 4.878787878787879,
                4.949494949494949, 5.02020202020202, 5.090909090909091, 5.161616161616162, 5.232323232323232,
                5.303030303030303, 5.373737373737374, 5.444444444444445, 5.515151515151515, 5.585858585858586,
                5.656565656565657, 5.727272727272727, 5.797979797979798, 5.8686868686868685, 5.9393939393939394,
                6.0101010101010095, 6.08080808080808, 6.151515151515151, 6.222222222222222, 6.292929292929292,
                6.363636363636363, 6.434343434343434, 6.505050505050505, 6.575757575757575, 6.646464646464646,
                6.717171717171717, 6.787878787878787, 6.858585858585858, 6.929292929292929, 7.0};

        TransferFunction tf6 = new TransferFunction(new double[]{1.0, 3.0, 3.0}, new double[]{1.0, 2.0, 1.0});
        double[] yOut = tf6.step(timePoints, new double[]{1.0, 0.0}).getResponse();

        System.out.println(Arrays.toString(yOut));

        TransferFunction tf7 = new TransferFunction(new double[]{1.0, 3.0, 3.0}, new double[]{1.0, 2.0, 1.0});
        double[] yOut2 = tf7.step(new double[]{0.0}).getResponse();

        System.out.println(Arrays.toString(yOut2));

        TransferFunction tf8 = new TransferFunction(new double[] {1, 1, 0}, new double[] {1, 8, 25, 12, 1, 9 , 8, 10});
        double[] wn = tf8.findFrequencies(9);

        System.out.println(Arrays.toString(wn));

        TransferFunction tf9 = new TransferFunction(new double[] {1, 0.1, 7.5}, new double[] {1, 0.12, 9, 0, 0});

        System.out.println(Arrays.toString(tf9.findFrequencies(9)));

        // Specs for band pass filter
        FilterSpecs.BandPassSpecs bpSpecs = new FilterSpecs.BandPassSpecs();
        // The bandwidth of the filter starts at the LowerPassBandFrequency and
        // ends at the UpperPassBandFrequency. The filter has lower stop band
        // which is set LowerStopBandFrequency and the upper stop band can be set
        // with UpperStopBandFrequency. The attenuation at the stop bands can be
        // set with the LowerStopBandAttenuation and UpperStopBandAttenuation
        // respectively. In a frequency spectrum, the order of the frequencies will be:
        // LowerStopBandFrequency < LowerPassBandFrequency < UpperPassBandFrequency <
        // UpperStopBandFrequency
        bpSpecs.setLowerPassBandFrequency(190.0); // 190 Hz lower pass band frequency
        bpSpecs.setUpperPassBandFrequency(210.0); // 210 Hz upper pass band frequency
        bpSpecs.setLowerStopBandFrequency(180.0); // 180 Hz lower stop band frequency
        bpSpecs.setUpperStopBandFrequency(220.0); // 220 Hz upper stop band frequency
        bpSpecs.setPassBandRipple(0.2); // 0.2 dB gain/ripple refer to note
        bpSpecs.setStopBandAttenuation(20.0); // 20 dB attenuation in the stop band

        FilterOrderResults.OrderAndCutoffFrequencies nW0W1 = Elliptic.ellipord(bpSpecs);
        TransferFunction el = Elliptic.newBandPass(nW0W1.getOrder(), bpSpecs.getPassBandRipple(),
                bpSpecs.getStopBandAttenuation(), nW0W1.getLowerCutoffFrequency(), nW0W1.getUpperCutoffFrequency());

        System.out.println(Arrays.toString(el.findFrequencies(9)));

        FilterSpecs.BandStopSpecs bsSpecs = new FilterSpecs.BandStopSpecs();
        // The notch of the filter starts at the LowerStopBandFrequency and
        // ends at the UpperStopBandFrequency. The filter has lower pass band
        // which is set LowerPassBandFrequency and the upper pass band can be set
        // with UpperPassBandFrequency. The attenuation at the notch can be
        // set with the StopBandAttenuation parameter and the attenuation/ripple
        // in the pass band can be set with the PassBandRipple parameter.
        // In a frequency spectrum, the order of the frequencies will be:
        // LowerPassBandFrequency < LowerStopBandFrequency < UpperStopBandFrequency <
        // UpperPassBandFrequency
        bsSpecs.setLowerPassBandFrequency(3.6e3); // 3600 Hz lower pass band frequency
        bsSpecs.setUpperPassBandFrequency(9.1e3); // 9100 Hz lower pass band frequency
        bsSpecs.setLowerStopBandFrequency(5.45e3); // 5450 Hz lower stop band frequency
        bsSpecs.setUpperStopBandFrequency(5.90e3); // 5900 Hz upper stop band frequency
        bsSpecs.setPassBandRipple(0.5); // 1.5 dB gain/ripple refer to note
        bsSpecs.setStopBandAttenuation(38.0); // 38 db attenuation at the notch

        nW0W1 = Elliptic.ellipord(bsSpecs);
        el = Elliptic.newBandStop(nW0W1.getOrder(), bsSpecs.getPassBandRipple(),
                bsSpecs.getStopBandAttenuation(), nW0W1.getLowerCutoffFrequency(), nW0W1.getUpperCutoffFrequency());

        System.out.println(Arrays.toString(el.findFrequencies(9)));

        TransferFunction tf10 = new TransferFunction(new double[] {1, 0.1, 7.5, 0, 1, 0, 0 ,-10},
                new double[] {1, 0, -20, 1e13, 0, 8, 1, 0.12, 9, 0, 0});

        System.out.println(Arrays.toString(tf10.findFrequencies(9)));
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
