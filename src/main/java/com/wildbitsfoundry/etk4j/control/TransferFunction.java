package com.wildbitsfoundry.etk4j.control;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.linearalgebra.EigenvalueDecomposition;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrices;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;
import com.wildbitsfoundry.etk4j.util.NumArrays;

/***
 *
 * @author StaticBeagle
 *
 */
public class TransferFunction {
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

    public double[] getPhaseAt(double[] frequencies) {
        double[] phase = new double[frequencies.length];
        for (int i = 0; i < frequencies.length; ++i) {
            phase[i] = this.evaluateAt(frequencies[i]).arg() * (180 / Math.PI);
        }
        unwrapPhase(phase);
        return phase;
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
        Complex[] mag = ComplexArrays.conv(a, ComplexArrays.conj(a));
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

        Complex[] conv = ComplexArrays.conv(num, conjDen);
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
            double phase = this.getPhaseAt(wcg);
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
    public double getPhaseAt(double f) {
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

    private TimeResponseResults lSim(double[][] input, double[] time, double[] initialConditions,
                                     StateSpace ss) {
        if(input.length != time.length) {
            throw new IllegalArgumentException("The input array and the time array must have the same length");
        }

        double[][] U = input;

        Matrix A = ss.getA();
        Matrix B = ss.getB();
        Matrix C = ss.getC();
        Matrix D = ss.getD();

        final int noStates = A.getRowCount();
        final int noInputs = B.getColumnCount();
        final int noSteps = time.length;

        // initial conditions
        double[] x0 = initialConditions == null ? new double[noStates] : initialConditions;
        double[][] xOut = new double[noSteps][noStates];

        if(time[0] == 0.0) {
            xOut[0] = x0;
        } else if(time[0] > 0.0) {
            xOut[0] = NumArrays.dot(x0, A.transpose().multiply(time[0]).expm().getAs2DArray());
        } else {
            throw new IllegalArgumentException("Initial time must be non negative");
        }

        // TODO
        if(noSteps == 1) {

        }

        double dt = time[1] - time[0];
        double[] delta = new double[time.length - 2];
        for(int i = 1; i < time.length - 1; ++i) {
            delta[i - 1] = (time[i + 1] - time[i]) / dt;
        }
        if(!NumArrays.allClose(delta, 1.0)) {
            throw new NonUniformTimeStepsException();
        }

        A.multiplyEquals(dt);
        B.multiplyEquals(dt);
        double[][] M = new double[noStates + noInputs][];
        for(int i = 0; i < M.length - 1; ++i) {
            M[i] = NumArrays.concatenate(A.getRow(i), B.getRow(i));
        }
        M[M.length - 1] = new double[noStates + noInputs];

        Matrix expMT = new Matrix(M).transpose().expm();
        double[][] Ad = new double[noStates][noStates];
        for(int i = 0; i < noStates; ++i) {
            double[] row = expMT.getRow(i);
            Ad[i] = Arrays.copyOf(row, row.length - 1);
        }
        double[][] Bd = expMT.subMatrix(noStates, noStates, 0, noStates - 1).getAs2DArray();
        for(int i = 1; i < noSteps; ++i) {
            xOut[i] = NumArrays.add(NumArrays.dot(xOut[i - 1], Ad), NumArrays.dot(U[i - 1], Bd));
        }
        double[] yOut = new double[noSteps];
        double[] c = C.transpose().getArray();
        double[] d = D.transpose().getArray();
        for(int i = 0; i < noSteps; ++i) {
            yOut[i] = NumArrays.dot(xOut[i], c) + NumArrays.dot(U[i], d);
        }
        return new TimeResponseResults(time, yOut, xOut);
    }

//    public StepResponse simulateTimeResponse(double[] input, double[] time) {
//        StateSpace ss = this.toStateSpace();
//        return lSim(input, time, null, ss);
//    }
//
//    public StepResponse simulateTimeResponse(double[] input, double[] time, double[] initialConditions) {
//        StateSpace ss = this.toStateSpace();
//        return lSim(input, time, initialConditions, ss);
//    }

    public StepResponse step() {
        return step(100);
    }

    public StepResponse step(int numberOfPoints) {
        return stepResponse(null, null, numberOfPoints);
    }

    public StepResponse step(double[] initialConditions, int numberOfPoints) {
        return stepResponse(null, initialConditions, numberOfPoints);
    }

    public StepResponse step(double[] time) {
        return stepResponse(time, null, time.length);
    }

    public StepResponse step(double[] time, double[] initialConditions) {
        return stepResponse(time, initialConditions, time.length);
    }

    private StepResponse stepResponse(double[] time, double[] initialConditions, int numberOfPoints) {
        StateSpace ss = this.toStateSpace();
        numberOfPoints = numberOfPoints == -1 ? 100 : numberOfPoints;
        time = time == null ? defaultResponseTimes(ss.getA(), numberOfPoints) : time;
        double[] input = NumArrays.ones(time.length);
        double[][] U = new double[time.length][1];
        for(int i = 0; i < U.length; ++i) {
            U[i][0] = input[i];
        }
        TimeResponseResults lSim = lSim(U, time, initialConditions, ss);
        return new StepResponse(lSim.getTime(), lSim.getResponse());
    }

    private double[] defaultResponseTimes(Matrix A, int numberOfPoints) {
        EigenvalueDecomposition eig = A.eig();
        double[] realEig = eig.getRealEigenvalues();
        for(int i = 0; i < realEig.length; ++i) {
            realEig[i] = Math.abs(realEig[i]);
        }
        double r = NumArrays.min(realEig);
        if(r == 0.0) {
            r = 1.0;
        }
        double tc = 1.0 / r;
        return NumArrays.linSpace(0.0, 7 * tc, numberOfPoints);
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

    public static void main(String[] args) {
        TransferFunction tf1 = new TransferFunction(new double[]{1000, 0}, new double[]{1, 25, 100, 9, 4});
        System.out.println(tf1);
        System.out.println(tf1.getMargins());

        double phase = tf1.getPhaseAt(70.4);
        System.out.println(phase);

        TransferFunction tf2 = new TransferFunction(new double[]{1, 3, 3}, new double[]{1, 2, 1});
        System.out.println(tf2.toStateSpace());

        TransferFunction tf3 = new TransferFunction(new double[]{5, 3, 4}, new double[]{8, 2, 9, 10});
        System.out.println(tf3.toStateSpace());

        TransferFunction tf4 = new TransferFunction(new double[]{5}, new double[]{3});
        System.out.println(tf4.toStateSpace());

        TransferFunction tf5 = new TransferFunction(new double[] {1.0, 3, 3}, new double[] {1.0, 2.0, 1});
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

        TransferFunction tf6 = new TransferFunction(new double[] {1.0, 3.0, 3.0}, new double[] {1.0, 2.0, 1.0});
        double[] yOut = tf6.step(timePoints, new double[] {1.0, 0.0}).getResponse();

        System.out.println(Arrays.toString(yOut));



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

    public void normalize() {
        this._rf.normalize();
    }
}
