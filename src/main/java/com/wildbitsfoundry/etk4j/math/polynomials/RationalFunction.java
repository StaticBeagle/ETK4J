package com.wildbitsfoundry.etk4j.math.polynomials;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.functions.ComplexUnivariateFunction;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;

public class RationalFunction implements UnivariateFunction, ComplexUnivariateFunction {
	private Polynomial _numerator;
	private Polynomial _denominator;

	public RationalFunction(Complex[] zeros, Complex[] poles) {
		this(zeros, poles, calculateGain(zeros, poles));
	}

	public RationalFunction(Complex[] zeros, Complex[] poles, double gain) {
		_numerator = new Polynomial(zeros);
		_denominator = new Polynomial(poles);

		_numerator.multiplyEquals(gain);
	}

	public RationalFunction(double num, Complex[] poles) {
		_numerator = new Polynomial(num);
		_denominator = new Polynomial(poles);

		_numerator.multiplyEquals(calculateGain(new Complex[] { Complex.fromReal(-1.0) }, poles));
	}

	public RationalFunction(RationalFunction rf) {
		_numerator = new Polynomial(rf._numerator);
		_denominator = new Polynomial(rf._denominator);
	}

	public RationalFunction(Polynomial num, Polynomial den) {
		_numerator = new Polynomial(num);
		_denominator = new Polynomial(den);
	}

	public RationalFunction(double num, Polynomial den) {
		_numerator = new Polynomial(num);
		_denominator = new Polynomial(den);
	}

	public RationalFunction(double[] num, double[] den) {
		_numerator = new Polynomial(num);
		_denominator = new Polynomial(den);
	}

	public Complex[] getZeros() {
		return _numerator.getRoots();
	}

	public Complex[] getPoles() {
		return _denominator.getRoots();
	}

	public Polynomial getNumerator() {
		return _numerator;
	}

	public Polynomial getDenominator() {
		return _denominator;
	}

	public RationalFunction add(final RationalFunction rf) {
		Polynomial numeratorLeftSide = _numerator.multiply(rf._denominator);
		Polynomial numeratorRightSide = _denominator.multiply(rf._numerator);
		Polynomial denominator = _denominator.multiply(rf._denominator);

		return new RationalFunction(numeratorLeftSide.add(numeratorRightSide), denominator);
	}

	public RationalFunction add(double scalar) {
		return new RationalFunction(_numerator.add(_denominator.multiply(scalar)), _denominator);
	}

	public RationalFunction subtract(RationalFunction rf) {
		Polynomial numeratorLeftSide = _numerator.multiply(rf._denominator);
		Polynomial numeratorRightSide = _denominator.multiply(rf._numerator);
		Polynomial denominator = _denominator.multiply(rf._denominator);

		return new RationalFunction(numeratorLeftSide.subtract(numeratorRightSide), denominator);
	}

	public RationalFunction multiply(final Polynomial p) {
		Polynomial numerator = _numerator.multiply(p);
		Polynomial denominator = new Polynomial(_denominator);
		return new RationalFunction(numerator, denominator);
	}

	private void multiplyEquals(Polynomial p) {
		_numerator.multiplyEquals(p);
	}

	public RationalFunction multiply(RationalFunction rf) {
		Polynomial numerator = _numerator.multiply(rf._numerator);
		Polynomial denominator = _denominator.multiply(rf._denominator);
		return new RationalFunction(numerator, denominator);
	}

	public RationalFunction multiply(double scalar) {
		Polynomial numerator = _numerator.multiply(scalar);
		Polynomial denominator = new Polynomial(_denominator);
		return new RationalFunction(numerator, denominator);
	}

	public RationalFunction substitute(double d) {
		Polynomial num = _numerator.substitute(d);
		Polynomial den = _denominator.substitute(d);

		return new RationalFunction(num, den);
	}

	public void substituteInPlace(double d) {
		_numerator.substituteInPlace(d);
		_denominator.substituteInPlace(d);
	}

	public RationalFunction substitute(Polynomial p) {
		Polynomial num = _numerator.substitute(p);
		Polynomial den = _denominator.substitute(p);

		return new RationalFunction(num, den);
	}

	public void substituteInPlace(Polynomial p) {
		_numerator.substituteInPlace(p);
		_denominator.substituteInPlace(p);
	}

	public RationalFunction substitute(RationalFunction rf) {
		RationalFunction result = new RationalFunction(this);
		subsOp(result, rf._numerator, rf._denominator);
		return result;
	}

	public void substituteInPlace(RationalFunction rf) {
		subsOp(this, rf._numerator, rf._denominator);
	}

	private static void subsOp(RationalFunction rf, Polynomial num, Polynomial den) {
		RationalFunction nump = rf._numerator.substitute(num, den);
		RationalFunction denp = rf._denominator.substitute(num, den);

		int numDegree = rf._numerator.degree();
		int denDegree = rf._denominator.degree();
		int diff = Math.abs(numDegree - denDegree);
		if (denDegree > numDegree) {
			nump.multiplyEquals(den.pow(diff));
		} else if (numDegree > denDegree) {
			denp.multiplyEquals(den.pow(diff));
		}
		rf._numerator = nump._numerator;
		rf._denominator = denp._numerator;
	}

	public double normalize() {
		double norm = _denominator.normalize();
		_numerator.multiplyEquals(1 / norm);
		return norm;
	}

	@Override
	public double evaluateAt(double x) {
		return _numerator.evaluateAt(x) / _denominator.evaluateAt(x);
	}

	@Override
	public Complex evaluateAt(Complex c) {
		return _numerator.evaluateAt(c).divide(_denominator.evaluateAt(c));
	}

	public Complex evaluateAt(double real, double imag) {
		Complex resultNum = _numerator.evaluateAt(real, imag);
		Complex resultDen = _denominator.evaluateAt(real, imag);
		resultNum.divideEquals(resultDen);
		return resultNum;
	}

	public static double calculateGain(Complex[] zeros, Complex[] poles) {
		// Compute gain k
		Complex knum = Complex.fromReal(1.0);
		for (Complex zero : zeros) {
			knum.multiplyEquals(zero.uminus());
		}
		Complex kden = Complex.fromReal(1.0);
		for (Complex pole : poles) {
			kden.multiplyEquals(pole.uminus());
		}
		kden.divideEquals(knum);
		return kden.real();
	}
}
