package com.wildbitsfoundry.etk4j.math.polynomials;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.wildbitsfoundry.etk4j.math.Formulas;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.functions.DifferentiableFunction;
import com.wildbitsfoundry.etk4j.math.functions.IntegrableFunction;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import com.wildbitsfoundry.etk4j.math.linearalgebra.EigenvalueDecomposition;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrices;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;
import com.wildbitsfoundry.etk4j.util.NumArrays;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkXYDimensions;

/**
 * 
 * @author Marcos L. Lopez-Rivera
 * @version
 *
 */
public class Polynomial implements UnivariateFunction, DifferentiableFunction, IntegrableFunction {

	protected double[] _coefs = null;
	protected Complex[] _roots = null;

	/***
	 * Constructs a polynomial with a constant value of 1.
	 */
	public Polynomial() {
		_coefs = new double[] { 1.0 };
	}

	/***
	 * Copy constructor
	 */
	public Polynomial(Polynomial polynomial) {
		int length = polynomial._coefs.length;
		this._coefs = new double[length];
		System.arraycopy(polynomial._coefs, 0, this._coefs, 0, length);
		if (polynomial._roots != null) {
			this._roots = ComplexArrays.deepCopy(polynomial._roots);
		}

	}

	/***
	 * Constructs a polynomial P(x) and initializes its coefficients to the
	 * coefficients passed as parameters. The coefficients are assumed to be in
	 * descending order i.e. [1, 3, 2] will generate P(x) = x^2 + 3x + 2
	 * 
	 * @param coefficients
	 *            Array of coefficients in descending order
	 */
	public Polynomial(double... coefficients) {

		int length = coefficients.length;
		int i = 0;
		while (i < length && coefficients[i] == 0.0) {
			++i;
		}
		this._coefs = i == length ? new double[] { 1.0 } : Arrays.copyOfRange(coefficients, i, length);

	}
	
	/***
	 * Constructs a polynomial P(x) and initializes its coefficients to the
	 * coefficients passed as parameters. The coefficients are assumed to be in
	 * descending order i.e. [1, 3, 2] will generate P(x) = x^2 + 3x + 2
	 * 
	 * @param coefficients
	 *            Array of coefficients in descending order
	 */
	public static Polynomial of(double... coefficients) {
		return new Polynomial(coefficients);
	}

	/***
	 * Creates a polynomial P(x) from an array of its roots. Say we have the
	 * following set of roots r = [ -1, -2 ] then: P(x) = (x + 1)*(x + 2) which
	 * is the solution to the polynomial: P(x) = x^2 + 3x + 2
	 * 
	 * @param roots
	 *            Array of roots
	 */
	public Polynomial(Complex... roots) {
		final int size = roots.length;

		Complex[] tmp = new Complex[size];
		Complex[] result = new Complex[size + 1];

		for (int i = 1; i <= size; ++i) {
			result[i] = new Complex();
		}
		result[0] = new Complex(1.0, 0.0);

		for (int i = 0; i < size; ++i) {
			// Fill up tmp
			for (int j = 0; j <= i; ++j) {
				tmp[j] = roots[i].multiply(result[j]);
			}
			for (int j = 0; j <= i; ++j) {
				result[j + 1].subtractEquals(tmp[j]);
			}
		}
		_coefs = new double[size + 1];
		for (int i = 0; i <= size; i++) {
			_coefs[i] = result[i].real();
		}
		_roots = ComplexArrays.deepCopy(roots);
	}
	
	/***
	 * 
	 * @return Returns the order of the polynomial
	 */
	public int degree() {
		return _coefs.length - 1;
	}

	/***
	 * Forces the lowest order coefficient of 
	 * the polynomial to be unity by dividing all the other coefficients by 
	 * the lowest order coefficient.
	 */
	public void denormalize() {
		int i;
		// Find the normalizing gain
		i = _coefs.length - 1;
		while (this._coefs[i] < Double.MIN_VALUE) {
			--i;
			if (i == 0)
				break;
		}
		for (int j = 0; j < this._coefs.length; j++) {
			this._coefs[j] = this._coefs[j] / this._coefs[i];
		}
	}

	/***
	 * Convert to monic Polynomial. Forces the highest order coefficient of 
	 * the polynomial to be unity by dividing all the other coefficients by 
	 * the highest order coefficient.
	 * 
	 * @return normalizing factor
	 */
	public double normalize() {
		double cn = this._coefs[0];
		for (int j = 1; j < this._coefs.length; j++) {
			this._coefs[j] = this._coefs[j] / cn;
		}
		this._coefs[0] = 1;
		return cn;
	}

	/***
	 * Multiply two polynomials
	 * 
	 * @param p
	 *            Another polynomial
	 * @return Pnew(x) = P(x) * poly
	 */
	public Polynomial multiply(final Polynomial p) {
		return new Polynomial(NumArrays.conv(_coefs, p._coefs));
	}
	
	public Polynomial multiply(double... coefs) {
		return new Polynomial(NumArrays.conv(_coefs, coefs));
	}

	/***
	 * Multiply two polynomials and stores the result
	 * 
	 * @param p
	 *            Another polynomial
	 */
	public void multiplyEquals(final Polynomial p) {
		_coefs = NumArrays.conv(_coefs, p._coefs);
		_roots = null;
	}
	
	/***
	 * Multiply two polynomials and stores the result
	 * 
	 * @param p
	 *            Another polynomial
	 */
	public void multiplyEquals(double... coefs) {
		_coefs = NumArrays.conv(_coefs, coefs);
		_roots = null;
	}

	public Polynomial multiply(double d) {
		return new Polynomial(NumArrays.multiply(_coefs, d));
	}

	public void multiplyEquals(double d) {
		for (int i = 0; i < this._coefs.length; i++) {
			this._coefs[i] = this._coefs[i] * d;
		}
	}

	public Polynomial add(final Polynomial p) {
		double[] result = addOp(this, p);
		return new Polynomial(result);
	}

	public void addEquals(final Polynomial p) {
		_coefs = addOp(this, p);
		_roots = null;
	}

	private static final double[] addOp(Polynomial p1, Polynomial p2) {
		final int p1Length = p1._coefs.length;
		final int p2Length = p2._coefs.length;
		double[] result = p1Length > p2Length ? Arrays.copyOf(p1._coefs, p1Length) : Arrays.copyOf(p2._coefs, p2Length);
		for (int i = p1Length - 1, j = p2Length - 1, k = result.length - 1; i >= 0 && j >= 0; --i, --j, --k) {
			result[k] = p1._coefs[i] + p2._coefs[j];
		}
		return result;
	}

	public Polynomial subtract(final Polynomial p) {
		double[] result = subtractOp(this, p);
		return new Polynomial(result);
	}

	public void subtractEquals(final Polynomial p) {
		_coefs = subtractOp(this, p);
		_roots = null;
	}

	private static final double[] subtractOp(Polynomial p1, Polynomial p2) {
		final int p1Length = p1._coefs.length;
		final int p2Length = p2._coefs.length;
		double[] result = null;
		if(p1Length > p2Length) {
			int diff = p1Length - p2Length;
			result = new double[p1Length];
			System.arraycopy(p2._coefs, 0, result, diff, p2Length);
			for (int i = 0; i < p1Length; ++i) {
				result[i] = p1._coefs[i] - result[i];
			}
		} else {
			int diff = p2Length - p1Length;
			result = new double[p2Length];
			System.arraycopy(p1._coefs, 0, result, diff, p1Length);
			for (int i = 0; i < p2Length; ++i) {
				result[i] = result[i] - p2._coefs[i];
			}
		}
		return result;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		final int length = this.degree();
		for (int i = 0; i <= length; ++i) {
			if (_coefs[i] == 0) {
				continue;
			}
			sb.append(_coefs[i] < 0 ? " - " : " + "); 
			int power = length - i;
			double coef = Math.abs(this._coefs[i]);
			if(i == length) {
				sb.append(String.format("%.4g", coef));
			}
			else if (coef != 1.0) {
				String x = power == 1 ? " * x" : " * x^" + power;
				sb.append(String.format("%.4g", coef)).append(x);
			} else {
				String x = power == 1 ? "x" : "x^" + power;
				sb.append(x);
			}
		}
		if(sb.length() > 3) {
			sb.replace(0, 3, _coefs[0] < 0 ? "-" : "");
		}
		return sb.toString();
	}

	public Polynomial derivative() {
		int length = _coefs.length - 1;
		double[] coefficients = new double[length];
		for (int i = 0; i < length; ++i) {
			coefficients[i] = _coefs[i] * (length - i);
		}
		return new Polynomial(coefficients);
	}

	/***
	 * Get polynomial coefficients
	 * 
	 * @return Array containing the coefficients in descending order
	 */
	public double[] getCoefficients() {
		return Arrays.copyOf(_coefs, _coefs.length);
	}

	/***
	 * Evaluates the polynomial at x using Horner's method
	 * 
	 * @param x
	 * @return The value of the polynomial at x
	 */
	@Override
	public double evaluateAt(double x) {
		// Horner's method
		double result = 0;
		for (double coefficient : this._coefs) {
			result = result * x + coefficient;
		}
		return result;
	}

	public void reverseInPlace() {
		_coefs = NumArrays.reverse(_coefs);
		_roots = null;
	}

	/***
	 * Evaluates the polynomial at real + j * imag using Horner's method
	 * 
	 * @param real
	 *            part of the complex number
	 * @param imaginary
	 *            part of the complex number
	 * @return The value of the polynomial at real + j * imag
	 */
	public Complex evaluateAt(double real, double imag) {
		// Horner's method
		Complex result = new Complex();
		for (double coef : _coefs) {
			result.multiplyEquals(real, imag);
			result.addEquals(coef);
		}
		return result;
	}
	
	public Complex evaluateAt(Complex c) {
		// Horner's method
		Complex result = new Complex();
		for (double coef : _coefs) {
			result.multiplyEquals(c);
			result.addEquals(coef);
		}
		return result;
	}

	public Complex[] getRoots() {
		// lazy creation of roots
		if (_roots == null) {
			int N = this.degree();
			_roots = new Complex[N];
			switch (N) {
			case 0:
				_roots = new Complex[0];
				break;
			case 1:
				_roots = new Complex[] { new Complex(-_coefs[1] / _coefs[0], 0) };
				break;
			case 2:
				_roots = Formulas.quadraticFormula(_coefs[0], _coefs[1], _coefs[2]);
				break;
			default:
				// Use generalized eigenvalue decomposition to find the roots
				Matrix c = Matrices.Companion(_coefs, N);
				EigenvalueDecomposition evd = c.eig();
				double[] realEig = evd.getRealEigenvalues();
				double[] imagEig = evd.getImagEigenvalues();
				for (int i = 0; i < N; i++) {
					_roots[i] = new Complex(realEig[i], imagEig[i]);
				}
			}
		}
		// Defensive copy
		return ComplexArrays.deepCopy(_roots);
	}

	/***
	 * Scales polynomial coefficients
	 * 
	 * @param d
	 */
	public void substituteInPlace(double d) {
		final int deg = this.degree();
		for (int i = 0; i < deg; ++i) {
			for (int j = i; j < deg; ++j) {
				_coefs[i] *= d;
			}
		}
		_roots = null;
	}

	/***
	 * Scales polynomial coefficients
	 * 
	 * @param d
	 */
	public Polynomial substitute(double d) {
		final int deg = this.degree();
		double[] result = Arrays.copyOf(_coefs, _coefs.length);
		for (int i = 0; i < deg; ++i) {
			for (int j = i; j < deg; ++j) {
				result[i] *= d;
			}
		}
		return new Polynomial(result);
	}

	public void substituteInPlace(Polynomial p) {
		Polynomial result = substituteOp(this, p);
		_coefs = result._coefs;
		_roots = null;
	}

	public Polynomial substitute(Polynomial p) {
		return substituteOp(this, p);
	}
	
	public RationalFunction substitute(final Polynomial num, final Polynomial den) {
		final int deg = this.degree();
		Polynomial nump = num.pow(deg);
		nump.multiplyEquals(_coefs[0]);
		
		// Pre-calculate powers 
		List<Polynomial> pows = new ArrayList<>(deg);
		pows.add(new Polynomial(1.0));
		pows.add(new Polynomial(num));
		for(int i = 2; i < deg; ++i) {
			pows.add(pows.get(i - 1).multiply(num));
		}
		
		Polynomial tmp = null;
		Polynomial denp = new Polynomial(1.0);
		for (int i = deg - 1, j = 1; i >= 0; --i, ++j) {
			tmp = pows.get(i);
			denp.multiplyEquals(den);
			tmp.multiplyEquals(denp);
			tmp.multiplyEquals(_coefs[j]);
			nump.addEquals(tmp);
		}
		return new RationalFunction(nump, denp);
	}

	private static Polynomial substituteOp(Polynomial src, Polynomial sub) {
		final int deg = src.degree();
		Polynomial result = sub.pow(deg);
		result.multiplyEquals(src._coefs[0]);
		Polynomial tmp = null;
		for (int i = deg - 1, j = 1; i >= 0; --i, ++j) {
			tmp = sub.pow(i);
			tmp.multiplyEquals(src._coefs[j]);
			result.addEquals(tmp);
		}
		return result;
	}

	public Polynomial pow(int n) {
		if (n < 0) {
			throw new IllegalArgumentException("Power must be >= 0");
		}
		if (n == 0) {
			return new Polynomial(new double[] { 1.0 });
		}
		double[] tmp = Arrays.copyOf(_coefs, _coefs.length);
		while (--n > 0) {
			tmp = NumArrays.conv(tmp, _coefs);
		}
		return new Polynomial(tmp);
	}

	public double getCoefficientAt(int index) {
		return _coefs[index];
	}
	
	/***
	 * Polynomial fit
	 * <pre>
	 * Finds a polynomial P(x) = c0 + c1*x + * c2*x^2 + ... + cn*x^n 
	 * of degree n that fits the data in y best in a least-square sense.
	 * </pre>
	 * @param x
	 *            Array of points for the independent variable x
	 * @param y
	 *            Array of solutions to y(x)
	 * @param n
	 *            Order of the polynomial
	 * @return Returns a polynomial of degree n fits the data y best in a
	 *         least-square sense
	 */
	public static Polynomial polyFit(double[] x, double[] y, int n) {
		checkXYDimensions(x, y);
		int dim = x.length;
		// Building the coefficient matrix
		Matrix A = Matrices.Vandermonde(x, dim, n + 1);
		// Building the solution vector
		Matrix b = new Matrix(y, dim);
		Matrix c = A.solve(b);

		double[] coeffs = new double[n + 1];
		for (int i = 0; i <= n; i++) {
			coeffs[i] = c.get(n - i, 0);
		}
		return Polynomial.of(coeffs);
	}

	public Polynomial integral() {
		return this.integral(0.0);
	}

	public Polynomial integral(double constant) {
		final int length = _coefs.length;
		double[] integral = new double[length + 1];
		for (int i = 0, j = length; i < length; ++i, --j) {
			integral[i] = (_coefs[i] / j);
		}
		integral[length] = constant;
		return new Polynomial(integral);
	}

	@Override
	public double integrate(double a, double b) {
		Polynomial integral = this.integral();
		return integral.evaluateAt(b) - integral.evaluateAt(a);
	}

	@Override
	public double differentiate(double x) {
		final int length = _coefs.length;
		double result = 0.0;
		for (int i = 0, j = 0; j < length; ++i, ++j) {
			result *= x;
			result += _coefs[i] * (length - j);
		}
		return result;
	}
}