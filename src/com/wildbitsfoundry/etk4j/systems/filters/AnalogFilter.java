package com.wildbitsfoundry.etk4j.systems.filters;

import com.wildbitsfoundry.etk4j.math.polynomials.Polynomial;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;
import com.wildbitsfoundry.etk4j.systems.TransferFunction;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.BandPassSpecs;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.LowPassSpecs;

public abstract class AnalogFilter {
	
	protected int _order;
	
	public static AnalogFilter newLowPassFilter(LowPassSpecs specs) {
		return specs._approxType.getLowPassPrototype().buildLowPass(specs);
	}
	
	public static AnalogFilter newHighPassFilter(HighPassSpecs specs) {
		return specs._approxType.getLowPassPrototype().toHighPass(specs);
	}
	
	public static AnalogFilter newBandPassFilter(BandPassSpecs specs) {
		return specs._approxType.getLowPassPrototype().toBandPass(specs);
	}
	
	public static AnalogFilter newBandStopFilter(BandStopSpecs specs) {
		return specs._approxType.getLowPassPrototype().toBandStop(specs);
	}
	
	public int getOrder() {
		return _order;
	}

	public static TransferFunction lpTobp(Polynomial num, Polynomial den, double w0, double bw) {
		Polynomial s = new Polynomial(1, 0);
		Polynomial s2w02 = new Polynomial(bw, 0, bw * w0 * w0);

		RationalFunction bp = new RationalFunction(num, den);
		bp.substituteInPlace(new RationalFunction(s2w02, s));

		bp.normalize();

		return new TransferFunction(bp);
	}

	public static TransferFunction lpTobs(Polynomial num, Polynomial den, double w0, double bw) {
		Polynomial s = new Polynomial(1, 0);
		Polynomial s2w02 = new Polynomial(bw, 0, bw * w0 * w0);

		RationalFunction bp = new RationalFunction(num, den);
		bp.substituteInPlace(new RationalFunction(s, s2w02));

		bp.normalize();

		return new TransferFunction(bp);
	}

	public static TransferFunction lpTohp(Polynomial numerator, Polynomial denominator) {
		return lpTohp(numerator.getCoefficients(), denominator.getCoefficients());
	}

	public static TransferFunction lpTohp(double[] numerator, double[] denominator) {
		final int numDegree = numerator.length - 1;
		final int denDegree = denominator.length - 1;
		final int filterOrder = numDegree + denDegree + 1;

		// Reverse coefficients then scale them by the
		// order of the denominator i.e. pad with zeros
		double[] hpNumerator = new double[filterOrder];
		for (int i = numDegree, j = 0; i >= 0; --i, ++j) {
			hpNumerator[j] = numerator[i];
		}

		// Reverse coefficients then scale them by the
		// order of the numerator i.e. pad with zeros
		double[] hpDenominator = new double[filterOrder];
		for (int i = denDegree, j = 0; i >= 0; --i, ++j) {
			hpDenominator[j] = denominator[i];
		}

		return new TransferFunction(hpNumerator, hpDenominator);
	}
}
