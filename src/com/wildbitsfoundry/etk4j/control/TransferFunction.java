package com.wildbitsfoundry.etk4j.control;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;

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

	public static void unwrapPhase(double[] phase) {
		int length = phase.length;
		double[] dp = new double[length];
		double[] dps = new double[length];
		double[] dpcorr = new double[length];
		double[] cumulativesum = new double[length];

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
			dpcorr[j] = dps[j] - dp[j];
		}
		// Ignore correction when incremental variation is smaller than cutoff
		for (j = 0; j < length - 1; j++) {
			if (Math.abs(dp[j]) < cutoff) {
				dpcorr[j] = 0;
			}
		}
		// Find cumulative sum of deltas
		cumulativesum[0] = dpcorr[0];
		for (j = 1; j < length - 1; j++) {
			cumulativesum[j] = cumulativesum[j - 1] + dpcorr[j];
		}
		// Integrate corrections and add to P to produce smoothed phase values
		for (j = 1; j < length; j++) {
			phase[j] += cumulativesum[j - 1];
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
			format = new String[] { " ", padding };
		} else if (numLength < denLength) {
			format = new String[] { padding, " " };
		} else {
			format = new String[] { " ", " " };
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
		
		for(int i = length - 1, j = 0; i >= 0; --i, ++j) {
			poly[i] *= getComplexSign(j);
			if(j % 2 == 0) {
				result[i] = new Complex(poly[i], 0.0);
			} else {
				result[i] = new Complex(0.0,poly[i]);
			}
		}
		return result;
	}
	
	private static int getComplexSign(int power) {
		int result = power - 4 * (power / 4);
		return result == 2 || result == 3 ? - 1 : 1;
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

		Complex[] solution = magPolyDen.subtract(magPolyNum).getRoots();
		double[] wgc = new double[solution.length];
		int j = 0;
		for (int i = 0; i < solution.length; i++) {
			double real = solution[i].real();
			if (solution[i].imag() == 0 &&  real > 0) {
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
		for(double f : freqs) {
			double slope = this.getMagnitudeAt(f) - this.getMagnitudeAt(f + 1e-12);
			if(slope > 0.0) {
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
		for(int i = 0; i < conv.length; ++i) {
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
		for(double f : freqs) {
			double slope = this.getMagnitudeAt(f) - this.getMagnitudeAt(f + 1e-12);
			if(slope > 0.0) {
				return f;
			}
		}
		return Double.POSITIVE_INFINITY;
	}

	public Margins getMargins() {
		double wcg = this.getGainCrossoverFrequency();
		double wcp = this.getPhaseCrossoverFrequency();

		double pm = wcg;
		if(wcg != Double.POSITIVE_INFINITY) {
			double phase = this.getPhaseAt(wcg);
			pm = 180.0 + phase;
		}
		double gm = wcp;
		if(wcp != Double.POSITIVE_INFINITY) {
			double mag = this.getMagnitudeAt(wcp);
			gm = 20 * Math.log10(1.0 / mag);
		}

		Margins m = new Margins();
		m.GainCrossOverFrequency = wcg;
		m.PhaseCrossOverFrequency = wcp;
		m.GainMargin = gm;
		m.PhaseMargin = pm;

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
		for(int i = 0; i < roots.length; ++i) {
			if(roots[i].imag() == 0) {
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

	public static void main(String[] args) {
		TransferFunction tf1 = new TransferFunction(new double[] { 1000, 0 }, new double[] { 1, 25, 100, 9, 4});
		System.out.println(tf1);
		System.out.println(tf1.getMargins());
		
		double phase = tf1.getPhaseAt(70.4);
		System.out.println(phase);
		
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
