package com.wildbitsfoundry.etk4j.math.polynomials;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;

public class RationalFunction implements UnivariateFunction {
	private Polynomial _numerator;
	private Polynomial _denominator;

	public RationalFunction(Complex[] zeros, Complex[] poles) {
		_numerator = new Polynomial(zeros);
		_denominator = new Polynomial(poles);
	}
	
	public RationalFunction(double num, Complex[] poles) {
		_numerator = new Polynomial(num);
		_denominator = new Polynomial(poles);
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
		return _numerator.roots();
	}

	public Complex[] getPoles() {
		return _denominator.roots();
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

	public RationalFunction multiply(RationalFunction rf) {
		Polynomial numerator = _numerator.multiply(rf._numerator);
		Polynomial denominator = _denominator.multiply(rf._denominator);
		return new RationalFunction(numerator, denominator);
	}

	public RationalFunction multiply(double scalar) {
		Polynomial numerator = _numerator.multiply(scalar);
		Polynomial denominator = _denominator.multiply(1.0);
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
		Polynomial num = subPolyOp(_numerator, rf);
		Polynomial den = subPolyOp(_denominator, rf);

		num.multiplyEquals(rf._denominator.pow(_numerator.degree()));
		den.multiplyEquals(rf._denominator.pow(_denominator.degree()));

		return new RationalFunction(num, den);
	}

	public void substituteInPlace(RationalFunction rf) {
		Polynomial num = subPolyOp(_numerator, rf);
		Polynomial den = subPolyOp(_denominator, rf);

		num.multiplyEquals(rf._denominator.pow(_denominator.degree()));
		den.multiplyEquals(rf._denominator.pow(_numerator.degree()));

		_numerator = num;
		_denominator = den;
	}
	
	public double normalize() {
		double norm = _denominator.normalize();
		_numerator.multiplyEquals(1 / norm);
		return norm;
	}

	private static Polynomial subPolyOp(Polynomial src, RationalFunction sub) {
		final int deg = src.degree();
		Polynomial result = sub._numerator.pow(deg);
		result.multiplyEquals(src._coefs[0]);
		Polynomial tmp = null;
		for (int i = deg - 1, j = 1; i >= 0; --i, ++j) {
			tmp = sub._numerator.pow(i);
			tmp.multiplyEquals(sub._denominator.pow(j));
			tmp.multiplyEquals(src._coefs[j]);
			result.addEquals(tmp);
		}
		return result;
	}

	@Override
	public double evaluateAt(double x) {
		return _numerator.evaluateAt(x) / _denominator.evaluateAt(x);
	}

	public Complex evaluateAt(double real, double imag) {
		Complex resultNum = _numerator.evaluateAt(real, imag);
		Complex resultDen = _denominator.evaluateAt(real, imag);
		resultNum.divideEquals(resultDen);
		return resultNum;

	}
}
