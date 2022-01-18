package com.wildbitsfoundry.etk4j.math.complex;

import com.wildbitsfoundry.etk4j.math.MathETK;

public class Complex implements Comparable<Complex> {

	private double _real;
	private double _imag;

	public Complex() {
		_real = 0.0;
		_imag = 0.0;
	}

	public Complex(double real, double imag) {
		_real = real;
		_imag = imag;
	}

	public Complex(Complex c) {
		_real = c._real;
		_imag = c._imag;
	}

	public static Complex newComplex(Complex c) {
		return new Complex(c._real, c._imag);
	}

	public static Complex fromReal(double d) {
		return new Complex(d, 0.0);
	}

	public static Complex fromImaginary(double d) {
		return new Complex(0.0, d);
	}

	/***
	 * Constructs a complex number from magnitude and phase angle
	 * 
	 * @param r
	 *            magnitude
	 * @param theta
	 *            phase angle in radians
	 * @return (r * cos(theta), r * sin(theta))
	 */
	public static Complex fromPolar(double r, double theta) {
		return new Complex(r * Math.cos(theta), r * Math.sin(theta));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int hash = 1;
		long temp;
		temp = Double.doubleToLongBits(_imag);
		hash = prime * hash + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(_real);
		hash = prime * hash + (int) (temp ^ (temp >>> 32));
		return hash;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (!(obj instanceof Complex))
			return false;
		Complex other = (Complex) obj;
		if (Double.doubleToLongBits(_imag) != Double.doubleToLongBits(other._imag))
			return false;
		if (Double.doubleToLongBits(_real) != Double.doubleToLongBits(other._real))
			return false;
		return true;
	}

	public boolean isReal() {
		return _imag == 0;
	}

	public boolean isClose(Complex other, double tol) {
		if(!MathETK.isClose(Double.doubleToLongBits(_imag), Double.doubleToLongBits(other._imag), tol)) {
			return false;
		}
		if(!MathETK.isClose(Double.doubleToLongBits(_real), Double.doubleToLongBits(other._real), tol)) {
			return false;
		}
		return true;
	}

	@Override
	public int compareTo(Complex obj) {
		final int SMALLER = -1;
		final int EQUAL = 0;
		final int GREATER = 1;
		if (this == obj) {
			return EQUAL;
		}
		if (this._real > obj._real) {
			return GREATER;
		}
		if (this._real < obj._real) {
			return SMALLER;
		}
		return Double.compare(this._imag, obj._imag);
	}

	public double real() {
		return _real;
	}

	public double imag() {
		return _imag;
	}

	public double abs() {
		return MathETK.hypot(_real, _imag);
	}

	/***
	 * 
	 * @return The angle in radians where the x-axis is in polar coordinates
	 */
	public double arg() {
		return Math.atan2(_imag, _real);
	}

	/***
	 * Norm
	 * 
	 * @return magnitude squared
	 */
	public double norm() {
		return _real * _real + _imag * _imag;
	}

	public Complex conj() {
		return new Complex(_real, -_imag);
	}

	public Complex invert() {
		Complex result = new Complex(_real, _imag);
		invertOp(result);
		return result;
	}

	public void invertEquals() {
		invertOp(this);
	}

	public Complex add(Complex c) {
		Complex result = new Complex(_real, _imag);
		addOp(result, c);
		return result;
	}

	public Complex add(double d) {
		Complex result = new Complex(_real, _imag);
		addOp(result, d);
		return result;
	}
	
	public Complex add(double real, double imag) {
		Complex result = new Complex(_real, _imag);
		addOp(result, real, imag);
		return result;
	}


	public void addEquals(Complex c) {
		addOp(this, c);
	}

	public void addEquals(double d) {
		addOp(this, d);
	}
	
	public void addEquals(double real, double imag) {
		addOp(this, real, imag);
	}

	public Complex subtract(Complex c) {
		Complex result = new Complex(_real, _imag);
		subtractOp(result, c);
		return result;
	}

	public Complex subtract(double d) {
		Complex result = new Complex(_real, _imag);
		subtractOp(result, d);
		return result;
	}

	public void subtractEquals(Complex c) {
		subtractOp(this, c);
	}

	public void subtractEquals(double d) {
		subtractOp(this, d);
	}

	public Complex multiply(Complex c) {
		Complex result = new Complex(_real, _imag);
		multiplyOp(result, c);
		return result;
	}

	public Complex multiply(double real, double imag) {
		Complex result = new Complex(_real, _imag);
		multiplyOp(result, real, imag);
		return result;
	}

	public Complex multiply(double d) {
		Complex result = new Complex(_real, _imag);
		multiplyOp(result, d);
		return result;
	}

	public void multiplyEquals(Complex c) {
		multiplyOp(this, c);
	}

	public void multiplyEquals(double real, double imag) {
		multiplyOp(this, real, imag);
	}

	public void multiplyEquals(double d) {
		multiplyOp(this, d);
	}

	public Complex divide(Complex c) {
		Complex result = c.invert();
		multiplyOp(result, this);
		return result;
	}

	public Complex divide(double d) {
		Complex result = new Complex(_real, _imag);
		multiplyOp(result, 1.0 / d);
		return result;
	}

	public void divideEquals(Complex c) {
		Complex result = c.invert();
		multiplyOp(this, result);
	}

	public void divideEquals(double d) {
		multiplyOp(this, 1.0 / d);
	}

	public Complex sqrt() {
		if (_real == 0 && _imag == 0) {
			return new Complex();
		}

		double z = Math.sqrt(0.5 * (Math.abs(_real) + this.abs()));
		if (_real >= 0) {
			return new Complex(z, 0.5 * _imag / z);
		} else {
			return new Complex(0.5 * Math.abs(_imag) / z, Math.copySign(z, _imag));
		}
	}

	public Complex pow(Complex c) {
		return this.log().multiply(c).exp();
	}

	public Complex pow(double d) {
		return this.log().multiply(d).exp();
	}

	public Complex pow2() {
		double real = _real * _real - _imag * _imag;
		double imag = 2 * _real * _imag;
		return new Complex(real, imag);
	}

	public Complex log() {
		return new Complex(Math.log(this.abs()), this.arg());
	}

	public Complex exp() {
		double exp = Math.exp(_real);
		return new Complex(exp * Math.cos(_imag), exp * Math.sin(_imag));
	}

	public Complex sqrt1z() {
		Complex result = Complex.fromReal(1.0);
		result.subtractEquals(this.pow2());
		return result.sqrt();
	}

	public Complex sin() {
		return new Complex(Math.sin(_real) * Math.cosh(_imag), Math.cos(_real) * Math.sinh(_imag));
	}

	public Complex asin() {
		return sqrt1z().add(this.multiply(0.0, 1.0)).log().multiply(0.0, -1.0);
	}

	public Complex cos() {
		return new Complex(Math.cos(_real) * Math.cosh(_imag), -Math.sin(_real) * Math.sinh(_imag));
	}

	public Complex acos() {
		return this.add(this.sqrt1z().multiply(0.0, 1.0)).log().multiply(0.0, -1.0);
	}

	public Complex tan() {
		if (_imag > 20.0) {
			return Complex.fromImaginary(1.0);
		}
		if (_imag < -20) {
			return Complex.fromImaginary(-1.0);
		}

		double dreal = 2.0 * _real;
		double dimag = 2.0 * _imag;

		double tmp = 1.0 / (Math.cos(dreal) + Math.cosh(dimag));
		return new Complex(Math.sin(dreal) * tmp, Math.sinh(dimag) * tmp);
	}

	public Complex atan() {
		Complex i = Complex.fromImaginary(1.0);
		return this.add(i).divide(i.subtract(this)).log().multiply(i.multiply(new Complex(0.5, 0.0)));
	}

	public Complex uminus() {
		return new Complex(-_real, -_imag);
	}

	@Override
	public String toString() {
		return String.format("(%.4f %s %.4fj)", _real, _imag >= 0.0 ? "+" : "-", Math.abs(_imag));
	}

	private static final void invertOp(Complex c) {
		double mag = 1.0 / c.norm();
		c._real *= mag;
		c._imag *= -mag;
	}

	private static final void addOp(Complex c, double d) {
		c._real += d;
	}
	
	private static final void addOp(Complex c, double real, double imag) {
		c._real += real;
		c._imag += imag;
	}

	private static final void addOp(Complex c1, Complex c2) {
		c1._real += c2._real;
		c1._imag += c2._imag;
	}

	private static final void subtractOp(Complex c, double d) {
		c._real -= d;
	}

	private static final void subtractOp(Complex c1, Complex c2) {
		c1._real -= c2._real;
		c1._imag -= c2._imag;
	}

	private static final void multiplyOp(Complex c1, Complex c2) {
		double re = c1._real * c2._real - c1._imag * c2._imag;
		c1._imag = c1._real * c2._imag + c1._imag * c2._real;
		c1._real = re;
	}

	private static final void multiplyOp(Complex c1, double real, double imag) {
		double re = c1._real * real - c1._imag * imag;
		c1._imag = c1._real * imag + c1._imag * real;
		c1._real = re;
	}

	private static final void multiplyOp(Complex c, double d) {
		c._real *= d;
		c._imag *= d;
	}

	public Complex sinh() {
		return new Complex(Math.sinh(_real) * Math.cos(_imag), Math.cosh(_real) * Math.sin(_imag));
	}

	public Complex cosh() {
		return new Complex(Math.cosh(_real) * Math.cos(_imag), Math.sinh(_real) * Math.sin(_imag));
	}

	public Complex tanh() {
		Complex num = new Complex(Math.tanh(_real), Math.tan(_imag));
		Complex den = new Complex(1.0, Math.tanh(_real) * Math.tan(_imag));
		return num.divide(den);
	}
}
