package com.wildbitsfoundry.etk4j.math.polynomials;

import com.wildbitsfoundry.etk4j.math.Formulas;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.functions.ComplexUnivariateFunction;
import com.wildbitsfoundry.etk4j.math.functions.DifferentiableFunction;
import com.wildbitsfoundry.etk4j.math.functions.IntegrableFunction;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import com.wildbitsfoundry.etk4j.math.linearalgebra.EigenvalueDecomposition;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;
import com.wildbitsfoundry.etk4j.util.DoubleArrays;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkXYDimensions;

/**
 * The {@code Polynomial} class represents a polynomial {@code P(x) = c<sub>0</sub>x<sup>n</sup> +
 * c<sub>1</sub>x<sup>n - 1</sup> + ... + c<sub>n</sub>x<sup>0</sup>}.
 */
public class Polynomial implements UnivariateFunction, ComplexUnivariateFunction, DifferentiableFunction,
        IntegrableFunction {

    protected double[] coefficients = null;
    protected Complex[] roots = null;

    /***
     * Copy constructor
     */
    public Polynomial(Polynomial polynomial) {
        int length = polynomial.coefficients.length;
        this.coefficients = new double[length];
        System.arraycopy(polynomial.coefficients, 0, this.coefficients, 0, length);
        if (polynomial.roots != null) {
            this.roots = ComplexArrays.deepCopy(polynomial.roots);
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
        while (coefficients[i] == 0.0) {
            if (++i == length) {
                break;
            }
        }
        this.coefficients = i == length ? new double[]{0.0} : Arrays.copyOfRange(coefficients, i, length);

    }

    /***
     * Creates a polynomial P(x) from an array of its roots. Say we have the
     * following set of roots r = [ -1, -2 ] then: P(x) = (x + 1)*(x + 2) which
     * is the solution to the polynomial: P(x) = x^2 + 3x + 2.
     *
     * @param roots Array of roots.
     */
    public Polynomial(Complex... roots) {
        List<Complex> finiteRoots = new ArrayList<>(roots.length);
        for (int i = 0; i < roots.length; ++i) {
            if (Math.abs(roots[i].real()) != Double.POSITIVE_INFINITY
                    && Math.abs(roots[i].imag()) != Double.POSITIVE_INFINITY) {
                finiteRoots.add(roots[i]);
            }
        }
        final int size = finiteRoots.size();

        Complex[] tmp = new Complex[size];
        Complex[] result = new Complex[size + 1];

        for (int i = 1; i <= size; ++i) {
            result[i] = new Complex();
        }
        result[0] = Complex.fromReal(1.0);

        for (int i = 0; i < size; ++i) {
            // Fill up tmp
            for (int j = 0; j <= i; ++j) {
                tmp[j] = finiteRoots.get(i).multiply(result[j]);
            }
            for (int j = 0; j <= i; ++j) {
                result[j + 1].subtractEquals(tmp[j]);
            }
        }
        coefficients = ComplexArrays.real(result);
        this.roots = ComplexArrays.deepCopy(roots);
    }

    /***
     *
     * @return Returns the order of the polynomial
     */
    public int degree() {
        return coefficients.length - 1;
    }

    /***
     * Forces the lowest order coefficient of
     * the polynomial to be unity by dividing all the other coefficients by
     * the lowest order coefficient.
     */
    public void denormalize() {
        int i;
        // Find the normalizing gain
        i = coefficients.length - 1;
        while (this.coefficients[i] == 0.0) {
            if (--i == 0)
                break;
        }
        for (int j = 0; j < this.coefficients.length; j++) {
            this.coefficients[j] = this.coefficients[j] / this.coefficients[i];
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
        double cn = 1.0 / this.coefficients[0];
        for (int j = 0; j < this.coefficients.length; j++) {
            this.coefficients[j] = this.coefficients[j] * cn;
        }
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
        return new Polynomial(DoubleArrays.convolution(coefficients, p.coefficients));
    }

    public Polynomial multiply(double... coefs) {
        return new Polynomial(DoubleArrays.convolution(coefficients, coefs));
    }

    /***
     * Multiply two polynomials and stores the result
     *
     * @param p
     *            Another polynomial
     */
    public void multiplyEquals(final Polynomial p) {
        coefficients = DoubleArrays.convolution(coefficients, p.coefficients);
        roots = null;
    }

    /***
     * Multiply two polynomials and stores the result
     *
     * @param coefs
     *            Another polynomial
     */
    public void multiplyEquals(double... coefs) {
        coefficients = DoubleArrays.convolution(coefficients, coefs);
        roots = null;
    }

    public Polynomial multiply(double d) {
        return new Polynomial(DoubleArrays.multiplyElementWise(coefficients, d));
    }

    public void multiplyEquals(double d) {
        for (int i = 0; i < this.coefficients.length; i++) {
            this.coefficients[i] = this.coefficients[i] * d;
        }
    }

    public Polynomial add(final Polynomial p) {
        double[] result = addOp(this, p);
        return new Polynomial(result);
    }

    public void addEquals(final Polynomial p) {
        coefficients = addOp(this, p);
        roots = null;
    }

    private static final double[] addOp(Polynomial p1, Polynomial p2) {
        final int p1Length = p1.coefficients.length;
        final int p2Length = p2.coefficients.length;
        double[] result = p1Length > p2Length ? Arrays.copyOf(p1.coefficients, p1Length) : Arrays.copyOf(p2.coefficients, p2Length);
        for (int i = p1Length - 1, j = p2Length - 1, k = result.length - 1; i >= 0 && j >= 0; --i, --j, --k) {
            result[k] = p1.coefficients[i] + p2.coefficients[j];
        }
        return result;
    }

    public Polynomial subtract(final Polynomial p) {
        double[] result = subtractOp(this, p);
        return new Polynomial(result);
    }

    public void subtractEquals(final Polynomial p) {
        coefficients = subtractOp(this, p);
        roots = null;
    }

    private static final double[] subtractOp(Polynomial p1, Polynomial p2) {
        final int p1Length = p1.coefficients.length;
        final int p2Length = p2.coefficients.length;
        double[] result = null;
        if (p1Length > p2Length) {
            int diff = p1Length - p2Length;
            result = new double[p1Length];
            System.arraycopy(p2.coefficients, 0, result, diff, p2Length);
            for (int i = 0; i < p1Length; ++i) {
                result[i] = p1.coefficients[i] - result[i];
            }
        } else {
            int diff = p2Length - p1Length;
            result = new double[p2Length];
            System.arraycopy(p1.coefficients, 0, result, diff, p1Length);
            for (int i = 0; i < p2Length; ++i) {
                result[i] = result[i] - p2.coefficients[i];
            }
        }
        return result;
    }

    @Override
    public String toString() {
        final int length = this.degree();
        if (length == 0 && coefficients[0] == 0.0) {
            return "0.000";
        }

        StringBuilder sb = new StringBuilder();
        for (int i = 0; i <= length; ++i) {
            if (coefficients[i] == 0) {
                continue;
            }
            sb.append(coefficients[i] < 0 ? " - " : " + ");
            int power = length - i;
            double coef = Math.abs(this.coefficients[i]);
            if (i == length) {
                sb.append(String.format("%.4g", coef));
            } else if (!MathETK.isClose(coef, 1.0, 1e-12)) {
                String x = power == 1 ? " * x" : " * x^" + power;
                sb.append(String.format("%.4g", coef)).append(x);
            } else {
                String x = power == 1 ? "x" : "x^" + power;
                sb.append(x);
            }
        }
        sb.replace(0, 3, coefficients[0] < 0 ? "-" : "");
        return sb.toString();
    }

    public Polynomial derivative() {
        int length = coefficients.length - 1;
        double[] coefficients = new double[length];
        for (int i = 0; i < length; ++i) {
            coefficients[i] = this.coefficients[i] * (length - i);
        }
        return new Polynomial(coefficients);
    }

    /***
     * Get polynomial coefficients
     *
     * @return Array containing the coefficients in descending order
     */
    public double[] getCoefficients() {
        return Arrays.copyOf(coefficients, coefficients.length);
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
        for (double coefficient : this.coefficients) {
            result = result * x + coefficient;
        }
        return result;
    }

    public void reverseInPlace() {
        coefficients = DoubleArrays.reverse(coefficients);
        roots = null;
    }

    /***
     * Evaluates the polynomial at real + j * imag using Horner's method
     *
     * @param real
     *            part of the complex number
     * @param imag
     *            part of the complex number
     * @return The value of the polynomial at real + j * imag
     */
    public Complex evaluateAt(double real, double imag) {
        // Horner's method
        Complex result = new Complex();
        for (double coef : coefficients) {
            result.multiplyEquals(real, imag);
            result.addEquals(coef);
        }
        return result;
    }

    @Override
    public Complex evaluateAt(Complex c) {
        // Horner's method
        Complex result = new Complex();
        for (double coefficient : coefficients) {
            result.multiplyEquals(c);
            result.addEquals(coefficient);
        }
        return result;
    }

    public Complex[] calculateRoots() {
        // lazy creation of roots
        if (roots == null) {
            int N = this.degree();
            switch (N) {
                case 0:
                    roots = new Complex[0];
                    break;
                case 1:
                    roots = new Complex[]{new Complex(-coefficients[1] / coefficients[0], 0)};
                    break;
                case 2:
                    roots = Formulas.quadraticFormula(coefficients[0], coefficients[1], coefficients[2]);
                    break;
                default:
                    // Use generalized eigenvalue decomposition to find the roots
                    roots = new Complex[N];
                    Matrix c = Matrix.companion(coefficients, N);
                    EigenvalueDecomposition evd = c.eig();
                    double[] realEig = evd.getRealEigenvalues();
                    double[] imagEig = evd.getImagEigenvalues();
                    roots = ComplexArrays.zip(realEig, imagEig);
            }
        }
        // Defensive copy
        return ComplexArrays.deepCopy(roots);
    }

    /***
     * Substitutes polynomial coefficients
     * For a polynomial P(x) = x^2 + x + 1, it substitutes the x by x * d
     * thus P(x * d) = x^2 * d^2 + x * d + 1
     *
     * @param d
     */
    public void substituteInPlace(double d) {
        roots = null;
        if (d == 0) {
            coefficients = new double[]{coefficients[this.degree()]};
            return;
        }

        final int deg = this.degree();
        for (int i = 0; i < deg; ++i) {
            for (int j = i; j < deg; ++j) {
                coefficients[i] *= d;
            }
        }
    }

    /***
     * Substitutes polynomial coefficients
     * For a polynomial P(x) = x^2 + x + 1, it substitutes the x by x * d
     * thus P(x * d) = x^2 * d^2 + x * d + 1
     * @param d
     * @return P(x * d)
     */
    public Polynomial substitute(double d) {
        final int deg = this.degree();
        double[] result = Arrays.copyOf(coefficients, coefficients.length);
        for (int i = 0; i < deg; ++i) {
            result[i] *= Math.pow(d, deg - i);
        }
        return new Polynomial(result);
    }

    /***
     * Substitutes polynomial coefficients with another polynomial
     * For a polynomial P(x) = 2 * x + 1, it substitutes the x by p(x).
     * Say p(x) = 3 * x + 2
     * then P(p(x)) = (3 * x + 2) * 2 + 1 = 6 * x + 5
     * @param p polynomial to be inserted into P(x)
     */
    public void substituteInPlace(Polynomial p) {
        Polynomial result = substituteOp(this, p);
        coefficients = result.coefficients;
        roots = null;
    }

    /***
     * Substitutes polynomial coefficients with another polynomial
     * For a polynomial P(x) = 2 * x + 1, it substitutes the x by p(x).
     * Say p(x) = 3 * x + 2
     * then P(p(x)) = (3 * x + 2) * 2 + 1 = 6 * x + 5
     * @param p polynomial to be inserted into P(x)
     * @return P(p ( x))
     */
    public Polynomial substitute(Polynomial p) {
        return substituteOp(this, p);
    }

    public RationalFunction substitute(final Polynomial num, final Polynomial den) {
        final int deg = this.degree();
        Polynomial nump = num.pow(deg);
        nump.multiplyEquals(coefficients[0]);

        // Pre-calculate powers
        List<Polynomial> pows = new ArrayList<>(deg);
        pows.add(new Polynomial(1.0));
        pows.add(new Polynomial(num));
        for (int i = 2; i < deg; ++i) {
            pows.add(pows.get(i - 1).multiply(num));
        }

        Polynomial tmp = null;
        Polynomial denp = new Polynomial(1.0);
        for (int i = deg - 1, j = 1; i >= 0; --i, ++j) {
            tmp = pows.get(i);
            denp.multiplyEquals(den);
            tmp.multiplyEquals(denp);
            tmp.multiplyEquals(coefficients[j]);
            nump.addEquals(tmp);
        }
        return new RationalFunction(nump, denp);
    }

    private static Polynomial substituteOp(Polynomial src, Polynomial sub) {
        final int deg = src.degree();
        Polynomial result = sub.pow(deg);
        result.multiplyEquals(src.coefficients[0]);
        Polynomial tmp = null;
        for (int i = deg - 1, j = 1; i >= 0; --i, ++j) {
            tmp = sub.pow(i);
            tmp.multiplyEquals(src.coefficients[j]);
            result.addEquals(tmp);
        }
        return result;
    }

    public Polynomial pow(int n) {
        if (n < 0) {
            throw new IllegalArgumentException("Power must be >= 0");
        }
        if (n == 0) {
            return new Polynomial(new double[]{1.0});
        }
        double[] tmp = Arrays.copyOf(coefficients, coefficients.length);
        while (--n > 0) {
            tmp = DoubleArrays.convolution(tmp, coefficients);
        }
        return new Polynomial(tmp);
    }

    public double getCoefficientAt(int index) {
        return coefficients[index];
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
        Matrix A = Matrix.vandermonde(x, dim, n + 1);
        // Building the solution vector
        Matrix b = new Matrix(y, dim);
        Matrix c = A.solve(b);

        double[] coeffs = new double[n + 1];
        for (int i = 0; i <= n; i++) {
            coeffs[i] = c.get(n - i, 0);
        }
        return new Polynomial(coeffs);
    }

    /**
     * Indefinite integral of the {@code Polynomial}.
     * This method is equivalent to calling {@link #integral()} with integration constant equal to 0.
     *
     * @return The indefinite integral of the polynomial.
     */
    public Polynomial integral() {
        return this.integral(0.0);
    }

    /**
     * Indefinite integral of the {@Polynomial}.
     *
     * @param constant The integration constant.
     * @return The indefinite integral of the polynomial.
     */
    public Polynomial integral(double constant) {
        final int length = coefficients.length;
        double[] integral = new double[length + 1];
        for (int i = 0, j = length; i < length; ++i, --j) {
            integral[i] = (coefficients[i] / j);
        }
        integral[length] = constant;
        return new Polynomial(integral);
    }

    /**
     * Definite integral of the {@code Polynomial}.
     *
     * @param a The lower bound of the integration.
     * @param b The upper bound of the integration.
     * @return The definite integral of the polynomial from a to b.
     */
    @Override
    public double integrate(double a, double b) {
        Polynomial integral = this.integral();
        return integral.evaluateAt(b) - integral.evaluateAt(a);
    }

    /**
     * Differentiate {@code Polynomial}.
     *
     * @param x The argument at which to evaluate the derivative of the polynomial.
     * @return The derivative of the polynomial evaluate at {@code x}.
     */
    @Override
    public double differentiate(double x) {
        final int length = coefficients.length - 1;
        double result = 0.0;
        for (int i = 0, j = 0; j < length; ++i, ++j) {
            result *= x;
            result += coefficients[i] * (length - j);
        }
        return result;
    }

    /**
     * Evaluate a {@Polynomial} from its coefficients.
     *
     * @param coefficients The coefficients of the polynomial to evaluate.
     * @param x            The argument at which to evaluate the polynomial.
     * @returnThe The value of the polynomial at {@code x}.
     */
    public static double polyval(double[] coefficients, double x) {
        return DoubleArrays.horner(coefficients, x);
    }

    /**
     * Evaluate a {@Polynomial} from its coefficients.
     *
     * @param coefficients The coefficients of the polynomial to evaluate.
     * @param x            The array of argument at which to evaluate the polynomial.
     * @returnThe The values of the polynomial at {@code x}.
     */
    public static double[] polyval(double[] coefficients, double[] x) {
        final int length = x.length;
        double[] result = new double[length];
        for (int i = 0; i < length; ++i) {
            result[i] = polyval(coefficients, x[i]);

        }
        return result;
    }

    /**
     * Evaluate a {@code Polynomial}.
     *
     * @param roots The roots of teh polynomial.
     * @param x     Argument at which to evaluate the polynomial.
     * @return The value of the polynomial at {@code x}.
     * @see <a href="https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.polyfromroots.html">
     * polyfromroots</a>
     */
    public static Complex polyvalFromRoots(Complex[] roots, Complex x) {
        return ComplexArrays.product(ComplexArrays.subtract(x, roots));
    }
}