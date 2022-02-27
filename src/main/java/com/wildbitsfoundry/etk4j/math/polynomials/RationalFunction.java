package com.wildbitsfoundry.etk4j.math.polynomials;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.functions.ComplexUnivariateFunction;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;
import com.wildbitsfoundry.etk4j.util.NumArrays;

public class RationalFunction implements UnivariateFunction, ComplexUnivariateFunction {
    private Polynomial numerator;
    private Polynomial denominator;

    public RationalFunction(Complex[] zeros, Complex[] poles) {
        this(zeros, poles, calculateZeroPoleGain(zeros, poles));
    }

    public RationalFunction(Complex[] zeros, Complex[] poles, double gain) {
        numerator = new Polynomial(zeros);
        denominator = new Polynomial(poles);

        numerator.multiplyEquals(gain);
    }

    public RationalFunction(RationalFunction rf) {
        numerator = new Polynomial(rf.numerator);
        denominator = new Polynomial(rf.denominator);
    }

    public RationalFunction(Polynomial num, Polynomial den) {
        numerator = new Polynomial(num);
        denominator = new Polynomial(den);
    }

    public RationalFunction(double[] num, double[] den) {
        double factor = 1.0 / den[0];
        num = NumArrays.multiplyElementWise(num, factor);
        den = NumArrays.multiplyElementWise(den, factor);
        numerator = new Polynomial(num);
        denominator = new Polynomial(den);
    }

    public Complex[] getZeros() {
        return numerator.calculateRoots();
    }

    public Complex[] getPoles() {
        return denominator.calculateRoots();
    }

    public Polynomial getNumerator() {
        return new Polynomial(numerator);
    }

    public Polynomial getDenominator() {
        return new Polynomial(denominator);
    }

    public RationalFunction add(final RationalFunction rf) {
        Polynomial numeratorLeftSide = numerator.multiply(rf.denominator);
        Polynomial numeratorRightSide = denominator.multiply(rf.numerator);
        Polynomial denominator = this.denominator.multiply(rf.denominator);

        return new RationalFunction(numeratorLeftSide.add(numeratorRightSide), denominator);
    }

    public RationalFunction add(double scalar) {
        return new RationalFunction(numerator.add(denominator.multiply(scalar)), denominator);
    }

    public RationalFunction subtract(RationalFunction rf) {
        Polynomial numeratorLeftSide = numerator.multiply(rf.denominator);
        Polynomial numeratorRightSide = denominator.multiply(rf.numerator);
        Polynomial denominator = this.denominator.multiply(rf.denominator);

        return new RationalFunction(numeratorLeftSide.subtract(numeratorRightSide), denominator);
    }

    public RationalFunction multiply(final Polynomial p) {
        Polynomial numerator = this.numerator.multiply(p);
        Polynomial denominator = new Polynomial(this.denominator);
        return new RationalFunction(numerator, denominator);
    }

    private void multiplyEquals(Polynomial p) {
        numerator.multiplyEquals(p);
    }

    public RationalFunction multiply(RationalFunction rf) {
        Polynomial numerator = this.numerator.multiply(rf.numerator);
        Polynomial denominator = this.denominator.multiply(rf.denominator);
        return new RationalFunction(numerator, denominator);
    }

    public RationalFunction multiply(double scalar) {
        Polynomial numerator = this.numerator.multiply(scalar);
        Polynomial denominator = new Polynomial(this.denominator);
        return new RationalFunction(numerator, denominator);
    }

    public RationalFunction substitute(double d) {
        Polynomial num = numerator.substitute(d);
        Polynomial den = denominator.substitute(d);

        return new RationalFunction(num, den);
    }

    public void substituteInPlace(double d) {
        numerator.substituteInPlace(d);
        denominator.substituteInPlace(d);
    }

    public RationalFunction substitute(Polynomial p) {
        Polynomial num = numerator.substitute(p);
        Polynomial den = denominator.substitute(p);

        return new RationalFunction(num, den);
    }

    public void substituteInPlace(Polynomial p) {
        numerator.substituteInPlace(p);
        denominator.substituteInPlace(p);
    }

    public RationalFunction substitute(RationalFunction rf) {
        RationalFunction result = new RationalFunction(this);
        subsOp(result, rf.numerator, rf.denominator);
        return result;
    }

    public void substituteInPlace(RationalFunction rf) {
        subsOp(this, rf.numerator, rf.denominator);
    }

    private static void subsOp(RationalFunction rf, Polynomial num, Polynomial den) {
        RationalFunction nump = rf.numerator.substitute(num, den);
        RationalFunction denp = rf.denominator.substitute(num, den);

        int numDegree = rf.numerator.degree();
        int denDegree = rf.denominator.degree();
        int diff = Math.abs(numDegree - denDegree);
        if (denDegree > numDegree) {
            nump.multiplyEquals(den.pow(diff));
        } else if (numDegree > denDegree) {
            denp.multiplyEquals(den.pow(diff));
        }
        rf.numerator = nump.numerator;
        rf.denominator = denp.numerator;
    }

    public double normalize() {
        double norm = denominator.normalize();
        numerator.multiplyEquals(norm);
        return norm;
    }

    @Override
    public double evaluateAt(double x) {
        return numerator.evaluateAt(x) / denominator.evaluateAt(x);
    }

    @Override
    public Complex evaluateAt(Complex c) {
        return numerator.evaluateAt(c).divide(denominator.evaluateAt(c));
    }

    public Complex evaluateAt(double real, double imag) {
        Complex resultNum = numerator.evaluateAt(real, imag);
        Complex resultDen = denominator.evaluateAt(real, imag);
        resultNum.divideEquals(resultDen);
        return resultNum;
    }

    private static double calculateZeroPoleGain(Complex[] zeros, Complex[] poles) {
        // Compute gain k
        Complex num = ComplexArrays.product(zeros).multiply(Math.pow(-1, zeros.length));
        Complex den = ComplexArrays.product(poles).multiply(Math.pow(-1, poles.length));
        den.divideEquals(num);
        return den.real();
    }

    public boolean isProper() {
        return numerator.degree() <= denominator.degree();
    }

    public boolean isStrictlyProper() {
        return numerator.degree() < denominator.degree();
    }
}
