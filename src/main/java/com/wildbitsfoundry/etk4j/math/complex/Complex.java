package com.wildbitsfoundry.etk4j.math.complex;

import com.wildbitsfoundry.etk4j.math.MathETK;

/**
 * The {@code Complex} class describes as complex number. <br>
 * The class internal storage is double for real and imaginary.
 */
public class Complex implements Comparable<Complex> {

    private double real;
    private double imag;

    /**
     * Constructs a {@code Complex} number with real and imaginary part equal to zero.
     */
    public Complex() {
        real = 0.0;
        imag = 0.0;
    }

    /**
     * Constructs a {@code Complex} number from the real and imaginary parts.
     *
     * @param real The real part of the complex number.
     * @param imag The imaginary part of the complex number.
     */
    public Complex(double real, double imag) {
        this.real = real;
        this.imag = imag;
    }

    /**
     * Copy Constructor
     * @param c The Complex number to be copied.
     */
    public Complex(Complex c) {
        real = c.real;
        imag = c.imag;
    }

    /**
     * Constructs a complex number with only real part and zero imaginary port.
     *
     * @param d The real part of the complex number.
     * @return A complex number with real part equal to {@code d} nd imaginary part equal to zero.
     */
    public static Complex fromReal(double d) {
        return new Complex(d, 0.0);
    }

    /**
     * Constructs a complex number with real part equal to zero and only imaginary port.
     *
     * @param d The imaginary part of the complex number.
     * @return A complex number with real part equal to zero and imaginary part equal to {@code d}.
     */
    public static Complex fromImaginary(double d) {
        return new Complex(0.0, d);
    }

    /**
     * Constructs a complex number from magnitude and phase angle
     *
     * @param r The magnitude of the complex number.
     * @param theta The phase angle in radians.
     * @return (r * cos ( theta), r * sin(theta))
     */
    public static Complex fromPolar(double r, double theta) {
        return new Complex(r * Math.cos(theta), r * Math.sin(theta));
    }

    /**
     * @see java.lang.Object#hashCode()
     */
    @Override
    public int hashCode() {
        final int prime = 31;
        int hash = 1;
        long temp;
        temp = Double.doubleToLongBits(imag);
        hash = prime * hash + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(real);
        hash = prime * hash + (int) (temp ^ (temp >>> 32));
        return hash;
    }

    /**
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
        if (Double.doubleToLongBits(imag) != Double.doubleToLongBits(other.imag))
            return false;
        if (Double.doubleToLongBits(real) != Double.doubleToLongBits(other.real))
            return false;
        return true;
    }

    /**
     * Is the complex number comprised of only real part.
     *
     * @return True if the imaginary part is equal to zero.
     */
    public boolean isReal() {
        return imag == 0;
    }

    /**
     * Compare magnitudes of complex numbers.
     * @see java.lang.Comparable(java.lang.Object)
     */
    @Override
    public int compareTo(Complex obj) {
        if (this == obj) {
            return 0;
        }
        return Double.compare(this.abs(), obj.abs());
    }

    /**
     * Real part of the complex number.
     *
     * @return The real part of the imaginary number.
     */
    public double real() {
        return real;
    }

    /**
     * Imaginary part of the complex number.
     *
     * @return The imaginary part of the complex number.
     */
    public double imag() {
        return imag;
    }

    /**
     * Absolute value of the complex number.
     *
     * @return The absolute value of the complex number. <br>
     * Theoretically this can be computed as {@code Math.sqrt(real * real + imag * imag)}. <br>
     * For further details, see {@link MathETK#hypot(double, double)}
     */
    public double abs() {
        return MathETK.hypot(real, imag);
    }

    /**
     * Argument of the complex number.
     *
     * @return The angle in radians where the x-axis is in polar coordinates.
     */
    public double arg() {
        return Math.atan2(imag, real);
    }

    /***
     * Norm of the complex number.
     *
     * @return The magnitude squared, {@code real * real + imag * imag}.
     */
    public double norm() {
        return real * real + imag * imag;
    }

    /**
     * Conjugate of the complex number.
     *
     * @return {@code new Complex(real, -imag)}.
     */
    public Complex conj() {
        return new Complex(real, -imag);
    }

    /**
     * Inverse of the complex number.
     *
     * @return {@code Complex a = 1 / new Complex(a)}.
     */
    public Complex invert() {
        Complex result = new Complex(real, imag);
        invertOp(result);
        return result;
    }

    /**
     * Inverse of the complex number in place.
     *
     * @return {@code Complex a = 1 / a. No new object is created.
     */
    public void invertEquals() {
        invertOp(this);
    }

    /**
     * Addition of complex numbers.
     *
     * @param c Summand.
     * @return The complex number + {@code c}.
     */
    public Complex add(Complex c) {
        Complex result = new Complex(real, imag);
        addOp(result, c);
        return result;
    }

    /**
     * Addition of a complex number and a real number.
     *
     * @param d Summand.
     * @return The complex number + {@code d}.
     */
    public Complex add(double d) {
        Complex result = new Complex(real, imag);
        addOp(result, d);
        return result;
    }

    /**
     * Addition of complex numbers.
     *
     * @param real Real part of the summand.
     * @param imag Imaginary part of the summand.
     * @return Complex number + {@code new Complex(real, imag}.
     */
    public Complex add(double real, double imag) {
        Complex result = new Complex(this.real, this.imag);
        addOp(result, real, imag);
        return result;
    }

    /**
     * Addition in place. <br>
     * Performs the equivalent of {@code Complex a += Complex c}.
     *
     * @param c Complex summand.
     */
    public void addEquals(Complex c) {
        addOp(this, c);
    }

    /**
     * Addition in place. <br>
     * Performs the equivalent of {@ code Complex a += d}.
     *
     * @param d Real summand.
     */
    public void addEquals(double d) {
        addOp(this, d);
    }

    /**
     * Addition in place. <br>
     * Performs the equivalent of {@code Complex a += new Complex(real, imag)}
     *
     * @param real the real part of the summand.
     * @param imag the imaginary part of the summand.
     */
    public void addEquals(double real, double imag) {
        addOp(this, real, imag);
    }

    /**
     * Subtraction of complex numbers.
     *
     * @param c Subtrahend.
     * @return The complex number - {@code c}.
     */
    public Complex subtract(Complex c) {
        Complex result = new Complex(real, imag);
        subtractOp(result, c);
        return result;
    }

    /**
     * Subtracion of a complex number and a real number.
     *
     * @param d Subtrahend.
     * @return The complex number - {@code d}.
     */
    public Complex subtract(double d) {
        Complex result = new Complex(real, imag);
        subtractOp(result, d);
        return result;
    }

    /**
     * Subtraction in place. <br>
     * Performs the equivalent of {@code Complex a -= Complex c}.
     *
     * @param c Complex subtrahend.
     */
    public void subtractEquals(Complex c) {
        subtractOp(this, c);
    }

    /**
     * Subtraction in place. <br>
     * Performs the equivalent of {@ code Complex a -= d}.
     *
     * @param d Real subtrahend.
     */
    public void subtractEquals(double d) {
        subtractOp(this, d);
    }

    public Complex multiply(Complex c) {
        Complex result = new Complex(real, imag);
        multiplyOp(result, c);
        return result;
    }

    public Complex multiply(double real, double imag) {
        Complex result = new Complex(this.real, this.imag);
        multiplyOp(result, real, imag);
        return result;
    }

    public Complex multiply(double d) {
        Complex result = new Complex(real, imag);
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
        Complex result = new Complex(real, imag);
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
        if (real == 0 && imag == 0) {
            return new Complex();
        }

        double z = Math.sqrt(0.5 * (Math.abs(real) + this.abs()));
        if (real >= 0) {
            return new Complex(z, 0.5 * imag / z);
        } else {
            return new Complex(0.5 * Math.abs(imag) / z, Math.copySign(z, imag));
        }
    }

    public Complex pow(Complex c) {
        return this.log().multiply(c).exp();
    }

    public Complex pow(double d) {
        return this.log().multiply(d).exp();
    }

    public Complex pow2() {
        double real = this.real * this.real - imag * imag;
        double imag = 2 * this.real * this.imag;
        return new Complex(real, imag);
    }

    public Complex log() {
        return new Complex(Math.log(this.abs()), this.arg());
    }

    public Complex exp() {
        double exp = Math.exp(real);
        return new Complex(exp * Math.cos(imag), exp * Math.sin(imag));
    }

    public Complex sqrt1z() {
        Complex result = Complex.fromReal(1.0);
        result.subtractEquals(this.pow2());
        return result.sqrt();
    }

    public Complex sin() {
        return new Complex(Math.sin(real) * Math.cosh(imag), Math.cos(real) * Math.sinh(imag));
    }

    public Complex asin() {
        return sqrt1z().add(this.multiply(0.0, 1.0)).log().multiply(0.0, -1.0);
    }

    public Complex cos() {
        return new Complex(Math.cos(real) * Math.cosh(imag), -Math.sin(real) * Math.sinh(imag));
    }

    public Complex acos() {
        return this.add(this.sqrt1z().multiply(0.0, 1.0)).log().multiply(0.0, -1.0);
    }

    public Complex tan() {
        if (imag > 20.0) {
            return Complex.fromImaginary(1.0);
        }
        if (imag < -20) {
            return Complex.fromImaginary(-1.0);
        }

        double dreal = 2.0 * real;
        double dimag = 2.0 * imag;

        double tmp = 1.0 / (Math.cos(dreal) + Math.cosh(dimag));
        return new Complex(Math.sin(dreal) * tmp, Math.sinh(dimag) * tmp);
    }

    public Complex atan() {
        Complex i = Complex.fromImaginary(1.0);
        return this.add(i).divide(i.subtract(this)).log().multiply(i.multiply(new Complex(0.5, 0.0)));
    }

    public Complex uminus() {
        return new Complex(-real, -imag);
    }

    @Override
    public String toString() {
        return String.format("(%.4f %s %.4fj)", real, imag >= 0.0 ? "+" : "-", Math.abs(imag));
    }

    private static final void invertOp(Complex c) {
        double mag = 1.0 / c.norm();
        c.real *= mag;
        c.imag *= -mag;
    }

    private static final void addOp(Complex c, double d) {
        c.real += d;
    }

    private static final void addOp(Complex c, double real, double imag) {
        c.real += real;
        c.imag += imag;
    }

    private static final void addOp(Complex c1, Complex c2) {
        c1.real += c2.real;
        c1.imag += c2.imag;
    }

    private static final void subtractOp(Complex c, double d) {
        c.real -= d;
    }

    private static final void subtractOp(Complex c1, Complex c2) {
        c1.real -= c2.real;
        c1.imag -= c2.imag;
    }

    private static final void multiplyOp(Complex c1, Complex c2) {
        double re = c1.real * c2.real - c1.imag * c2.imag;
        c1.imag = c1.real * c2.imag + c1.imag * c2.real;
        c1.real = re;
    }

    private static final void multiplyOp(Complex c1, double real, double imag) {
        double re = c1.real * real - c1.imag * imag;
        c1.imag = c1.real * imag + c1.imag * real;
        c1.real = re;
    }

    private static final void multiplyOp(Complex c, double d) {
        c.real *= d;
        c.imag *= d;
    }

    public Complex sinh() {
        return new Complex(Math.sinh(real) * Math.cos(imag), Math.cosh(real) * Math.sin(imag));
    }

    public Complex cosh() {
        return new Complex(Math.cosh(real) * Math.cos(imag), Math.sinh(real) * Math.sin(imag));
    }

    public Complex tanh() {
        Complex num = new Complex(Math.tanh(real), Math.tan(imag));
        Complex den = new Complex(1.0, Math.tanh(real) * Math.tan(imag));
        return num.divide(den);
    }
}
