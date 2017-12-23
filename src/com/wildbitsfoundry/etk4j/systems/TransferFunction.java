package com.wildbitsfoundry.etk4j.systems;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;
import com.wildbitsfoundry.etk4j.util.NumArrays;

public class TransferFunction {
	private RationalFunction _rf;

	private static final double radToDeg = 57.295779513082320;
	//private static final double rad2Hz = 6.283185307179586;

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

	public double getMagnitudeAt(double f) {
		return _rf.evaluateAt(0.0, f).abs();
	}

	public double[] getMagnitudeAt(final double[] f) {
		int length = f.length;
		double[] magnitude = new double[length];
		for (int i = 0; i < length; ++i) {
			magnitude[i] = this.getMagnitudeAt(f[i]);
		}
		return magnitude;
	}

	public double getPhaseAt(double f) {
		return _rf.evaluateAt(0.0, f).arg() * radToDeg;
	}
	
	public double[] getPhaseWrappedAt(double[] f) {
		int length = f.length;
		double[] phase = new double[length];

		for (int i = 0; i < length; ++i) {
			phase[i] = this.getPhaseAt(f[i]);
		}
		return phase;
	}

	public double[] getPhaseAt(final double[] f) {
		double[] phase = this.getPhaseWrappedAt(f);
		unwrapPhase(phase);
		return phase;
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
		Polynomial numerator = _rf.getNumerator();
		Polynomial denominator = _rf.getDenominator();
		
		int lengthNumerator = numerator.toString().length();
		int lengthDenominator = denominator.toString().length();
		java.lang.StringBuilder sb = new java.lang.StringBuilder();
		char[] divider = new char[Math.max(lengthNumerator, lengthDenominator) + 2];
		Arrays.fill(divider, '-');

		char[] padding = new char[(int) Math
				.floor((divider.length - Math.min(lengthNumerator, lengthDenominator)) * 0.5)];
		Arrays.fill(padding, ' ');
		if (lengthNumerator > lengthDenominator) {
			sb.append(" ").append(numerator.toString()).append(System.getProperty("line.separator"))
					.append(divider).append(System.getProperty("line.separator")).append(padding)
					.append(denominator.toString());
			return sb.toString().replace('x', 's');
		} else if (lengthNumerator < lengthDenominator) {
			sb.append(padding).append(numerator.toString()).append(System.getProperty("line.separator"))
					.append(divider).append(System.getProperty("line.separator")).append(" ")
					.append(denominator.toString());
			return sb.toString().replace('x', 's');
		} else {
			sb.append(" ").append(numerator.toString()).append(System.getProperty("line.separator"))
					.append(divider).append(System.getProperty("line.separator")).append(" ")
					.append(denominator.toString());
			return sb.toString().replace('x', 's');
		}
	}

	private static Polynomial getPolynomialMagnitude(Polynomial polynomial) {
		double[] num = polynomial.getCoefficients();
		Complex[] numjw = evalAtjw(num);
		return new Polynomial(norm(numjw));
	}
	
	private static double[] norm(Complex[] vector) {
		Complex[] mag = conv(vector, conj(vector));
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
	
	private static Complex[] conj(Complex[] poly) {
		Complex[] result = new Complex[poly.length];
		for(int i = 0; i < poly.length; ++i) {
			result[i] = poly[i].conj();
		}
		return result;
	}
	
	private static Complex[] conv(Complex[] a, Complex[] b) 
	{
		Complex[] result = new Complex[a.length + b.length - 1];
	
	       for (int i = 0; i < result.length; ++i) 
	       {
	    	   result[i] = new Complex();
	           for (int j = Math.max(0, i + 1 - b.length); j < Math.min(a.length, i + 1); ++j)
	           {
	        	   result[i] = result[i].add(a[j].multiply(b[i-j]));
	           }
	       }
	
	       return result;
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

		Complex[] solution = magPolyDen.subtract(magPolyNum).roots();
		double[] wgc = new double[solution.length];
		int j = 0;
		for (int i = 0; i < solution.length; i++) {
			if (Double.compare(solution[i].imag(), 0.0) == 0 && solution[i].real() > 0) {
				wgc[j] = solution[i].real();
				j++;
			}
		}
		return j > 0 ? Arrays.copyOfRange(wgc, 0, j) : new double[0];
	}

	public double getGainCrossoverFrequency() {
		double[] result = this.getAllGainCrossoverFrequencies();
		return result.length > 0 ? result[0] : Double.NaN;
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
		Complex[] conjDen = conj(evalAtjw(_rf.getDenominator().getCoefficients()));
		
		Complex[] conv = conv(num, conjDen);
		double[] imag = new double[conv.length];
		for(int i = 0; i < conv.length; ++i) {
			imag[i] = conv[i].imag();
		}
		
		Complex[] solution = new Polynomial(imag).roots();

		double[] wpc = new double[solution.length];
		int j = 0;
		for (int i = 0; i < solution.length; i++) {
			if (solution[i].imag() == 0 && solution[i].real() > 0) {
				wpc[j] = solution[i].real();
				j++;
			}
		}
		return j > 0 ? Arrays.copyOfRange(wpc, 0, j) : new double[0];
	}

	public double getPhaseCrossoverFrequency() {
		double[] result = this.getAllPhaseCrossoverFrequencies();
		return result.length > 0 ? result[0] : Double.NaN;
	}

	public Margins getMargins() {
		double[] wcg = this.getAllGainCrossoverFrequencies();
		double[] wcp = this.getAllPhaseCrossoverFrequencies();

		Arrays.sort(wcg);
		Arrays.sort(wcp);

		double pm = wcg[0] == Double.NaN ? Double.NaN : 180 + this.getPhaseAt(wcg[0]);
		double gm = wcp[0] == Double.NaN ? Double.NaN : 20 * Math.log10(1.0 / this.getMagnitudeAt(wcp[0]));

		return new Margins(wcg, wcp, pm, gm);
	}
	
	public void substituteInPlace(double d) {
		_rf.substituteInPlace(d);
	}

	public static void main(String[] args) {
		double[] logspace = NumArrays.logspace(-1, 1, 100);

		TransferFunction tf1 = new TransferFunction(new double[] { 1, 10, 1000 }, new double[] { 1, 25, 100, 9, 4 });
		System.out.println(tf1);
		System.out.println(tf1.getMargins());
		
		double[] phase = tf1.getPhaseAt(logspace);
		System.out.println(Arrays.toString(phase));
		for(int i = 0; i < phase.length; ++i) {
			System.out.printf("%.4f %.4f%n", logspace[i], phase[i]);
		}
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
