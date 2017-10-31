package com.wildbitsfoundry.etk4j.math.polynomials;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.interpolation2d.Spline2d;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;
import com.wildbitsfoundry.etk4j.util.ArrayUtils;

public class Polynomial2d {
	private double[] _coefs;
	private int _n;
	private int _m;

	public Polynomial2d(double[] coefs, int n, int m) {
		_coefs = Arrays.copyOf(coefs, coefs.length);
		_n = n;
		_m = m;
	}

	/***
	 * Finds the coefficients of a polynomial P(x,y) of degree n + m that fits
	 * the data z best in a least-square sense. Always make sure that the number
	 * of points for x,y,z are the same and that the total number of points >=
	 * (n + 1) * (m + 1)
	 * 
	 * @param x
	 *            Array of points for the independent variable x
	 * @param y
	 *            Array of points for the independent variable y
	 * @param z
	 *            Array of solutions to z(x,y)
	 * @param n
	 *            Maximum order allowed for x
	 * @param m
	 *            Maximum order allowed for y
	 * @return Returns a bivariate polynomial of degree n + m fits the data z
	 *         best in a least-square sense
	 */
	public static Polynomial2d polyFit2D(double[] x, double[] y, double[] z, int n, int m) {
		double[][] vars = new double[2][];
		vars[0] = x;
		vars[1] = y;

		int[] nm = new int[] { n, m };
		double[] coefs = polyfitN(vars, z, nm);
		return new Polynomial2d(coefs, n, m);
	}

	/***
	 * Finds the coefficients of a polynomial P(x,y) of degree n + m that fits
	 * the data z best in a least-square sense. Always make sure that the number
	 * of points for x,y,z are the same and that the total number of points >=
	 * (n + 1) * (m + 1) in a regular grid
	 * 
	 * @param x
	 *            Array of points for the independent variable x
	 * @param y
	 *            Array of points for the independent variable y
	 * @param z
	 *            Array of solutions to z(x,y)
	 * @param n
	 *            Maximum order allowed for x
	 * @param m
	 *            Maximum order allowed for y
	 * @return Returns a bivariate polynomial of degree n + m fits the data z
	 *         best in a least-square sense
	 */
	public static Polynomial2d polyFit2D(double[] x, double[] y, double[][] z, int n, int m) {
		double[][] vars = new double[2][];
		// Build grid
		x = ArrayUtils.repeat(x, z[0].length);

		y = ArrayUtils.repeat(y, z.length);
		Arrays.sort(y);
		vars[0] = x;
		vars[1] = y;

		int[] nm = new int[] { n, m };
		double[] coefs = polyfitN(vars, ArrayUtils.flatten(z), nm);
		return new Polynomial2d(coefs, n, m);
	}

	private static double[] polyfitN(double[][] vars, double[] z, int[] n) {

		// int order = (m + 1) * (n + 1); // ... * (nx + 1)
		int[] degrees = n;
		int dims = vars.length;
		// double[][] vars = new double[dims][];
		// vars[0] = x;
		// vars[1] = y;
		// ... vars[n]

		int k = dims - 1;
		double[][] vandermonde = new double[z.length][];
		for (int i = 0; i < z.length; ++i) {
			vandermonde[i] = buildPowerSeriesDescending(vars[k][i], degrees[k]);
			for (int j = k - 1; j >= 0; --j) {
				double[] tmp = buildPowerSeriesDescending(vars[j][i], degrees[j]);
				vandermonde[i] = ArrayUtils.kron(vandermonde[i], tmp);
			}
		}

		Matrix A = new Matrix(vandermonde);
		Matrix b = new Matrix(z, z.length);

		Matrix c = A.solve(b);

		return c.getArrayCopy();
	}

	private static double[] buildPowerSeriesDescending(double x, int power) {
		double[] series = new double[power + 1];
		int i = 0;
		while (power > 0) {
			series[i++] = Math.pow(x, power--);
		}
		series[i] = 1.0;
		return series;
	}

	public double[] evaluateAt(double[] x, double[] y) {
		if (x.length != y.length) {
			throw new IllegalArgumentException("The length of x and y arrays must the same");
		}
		double[] result = new double[x.length];
		final int n = _n;
		final int m = _m;

		for (int i = 0; i < x.length; ++i) {
			result[i] = _coefs[0];
		}

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < x.length; ++j) {
				result[j] = result[j] * x[j] + _coefs[i + 1];
			}
		}

		for (int i = 0; i < m; ++i) {
			int mj = (n + 1) * (i + 1);
			double[] g = new double[x.length];

			for (int j = 0; j < x.length; ++j) {
				g[j] = _coefs[mj];
			}

			for (int k = 0; k < n; ++k) {
				for (int j = 0; j < x.length; ++j) {
					g[j] = g[j] * x[j] + _coefs[mj + k + 1];
				}
			}

			for (int l = 0; l < x.length; l++) {
				result[l] = result[l] * y[l] + g[l];
			}
		}
		return result;
	}

	/***
	 * Evaluate a 3-D polynomial using Horner's method.
	 * 
	 * @param poly
	 *            Polynomial to be evaluated at the points polyecified by x and y,
	 *            which must have the same dimensions.
	 * @param x
	 *            Arrays of points for the independent variable x
	 * @param y
	 *            Arrays of points for the independent variable y
	 * @param n
	 *            Order of x
	 * @param m
	 *            order of y
	 * @return Array of values resulting from evaluating poly(x,y) at the points
	 *         polyecified by x and y
	 */
	public double evaluateAt(double x, double y) {
		double result = _coefs[0];
		final int n = _n;
		final int m = _m;

		for (int i = 0; i < n; ++i) {
			result = result * x + _coefs[i + 1];
		}

		for (int i = 0; i < m; ++i) {
			int mj = (n + 1) * (i + 1);
			double g = _coefs[mj];

			for (int k = 0; k < n; ++k) {
				g = g * x + _coefs[mj + k + 1];
			}

			result = result * y + g;
		}
		return result;
	}

	@Override
	public String toString() {
		java.lang.StringBuilder sb = new java.lang.StringBuilder();

		for (int i = this._m, k = 0; i >= 0; i--) {
			for (int j = this._n; j >= 0; j--) {
				if (j == 0 && i == 0) {
					sb.append(String.format("%.4e", _coefs[k]));
				} else if (j == 0) {
					if (i == 1) {
						sb.append(String.format("%.4e*Y + ", _coefs[k], i));
					} else {
						sb.append(String.format("%.4e*Y^%d + ", _coefs[k], i));
					}
				} else if (i == 0) {
					if (j == 1) {
						sb.append(String.format("%.4e*X + ", _coefs[k]));
					} else {
						sb.append(String.format("%.4e*X^%d + ", _coefs[k], j));
					}
				} else {
					if (i == 1 && j == 1) {
						sb.append(String.format("%.4e*X*Y + ", _coefs[k]));
					} else if (j == 1) {
						sb.append(String.format("%.4e*X*Y^%d + ", _coefs[k], i));
					} else if (i == 1) {
						sb.append(String.format("%.4e*X^%d*Y + ", _coefs[k], j));
					} else {
						sb.append(String.format("%.4e*X^%d*Y^%d + ", _coefs[k], j, i));
					}

				}
				k++;
			}
		}
		return sb.toString();
	}

	public static void main(String[] args) {
		double[] x = new double[] { 0, 1, 1, 2, 3, 4, 5, 6, 8, 9 };
		double[] y = new double[] { 1, 2, 7, 9, 10, 14, 15, 33, 48, 59 };

		double[] z = new double[] { 1, 1, 6, 3, 4, 9, 6, 8, 9 };

		System.out.println(Arrays.toString(polyFit2D(x, y, z, 2, 2)._coefs));

		double[] xp = new double[] { 1, 2, 3, 4, 5, 6, 7, 8 };
		double[] yp = new double[] { 1, 2, 3, 4, 5, 6, 7, 8 };

		double[][] zp = new double[][] { { 1, 4, 9, 16, 25, 36, 49, 64 }, { 4, 16, 36, 64, 100, 144, 196, 256 },
				{ 9, 36, 81, 144, 225, 324, 441, 576 }, { 16, 64, 144, 256, 400, 576, 784, 1024 },
				{ 25, 100, 225, 400, 625, 900, 1225, 1600 }, { 36, 144, 324, 576, 900, 1296, 1764, 2304 },
				{ 49, 196, 441, 784, 1225, 1764, 2401, 3136 }, { 64, 256, 576, 1024, 1600, 2304, 3136, 4096 } };

		
		
		Polynomial2d poly = polyFit2D(xp, yp, zp, 2, 2);
		System.out.println(Arrays.toString(poly._coefs));
		
		System.out.println(poly.evaluateAt(1, 1));
		System.out.println(poly.evaluateAt(2, 1));
		System.out.println(poly.evaluateAt(3, 1));
		System.out.println(poly.evaluateAt(4, 1));
		System.out.println(poly.evaluateAt(5, 1));
		System.out.println(poly.evaluateAt(6, 1));
		System.out.println(poly.evaluateAt(7, 1));
		System.out.println(poly.evaluateAt(8, 1));

		System.out.println(poly.evaluateAt(1, 2));
		System.out.println(poly.evaluateAt(2, 2));
		System.out.println(poly.evaluateAt(3, 2));
		System.out.println(poly.evaluateAt(4, 2));
		System.out.println(poly.evaluateAt(5, 2));
		System.out.println(poly.evaluateAt(6, 2));
		System.out.println(poly.evaluateAt(7, 2));
		System.out.println(poly.evaluateAt(8, 2));

		System.out.println(poly.evaluateAt(1, 3));
		System.out.println(poly.evaluateAt(2, 3));
		System.out.println(poly.evaluateAt(3, 3));
		System.out.println(poly.evaluateAt(4, 3));
		System.out.println(poly.evaluateAt(5, 3));
		System.out.println(poly.evaluateAt(6, 3));
		System.out.println(poly.evaluateAt(7, 3));
		System.out.println(poly.evaluateAt(8, 3));

		System.out.println(poly.evaluateAt(1, 4));
		System.out.println(poly.evaluateAt(2, 4));
		System.out.println(poly.evaluateAt(3, 4));
		System.out.println(poly.evaluateAt(4, 4));
		System.out.println(poly.evaluateAt(5, 4));
		System.out.println(poly.evaluateAt(6, 4));
		System.out.println(poly.evaluateAt(7, 4));
		System.out.println(poly.evaluateAt(8, 4));

		System.out.println(poly.evaluateAt(1, 5));
		System.out.println(poly.evaluateAt(2, 5));
		System.out.println(poly.evaluateAt(3, 5));
		System.out.println(poly.evaluateAt(4, 5));
		System.out.println(poly.evaluateAt(5, 5));
		System.out.println(poly.evaluateAt(6, 5));
		System.out.println(poly.evaluateAt(7, 5));
		System.out.println(poly.evaluateAt(8, 5));

		System.out.println(poly.evaluateAt(1, 6));
		System.out.println(poly.evaluateAt(2, 6));
		System.out.println(poly.evaluateAt(3, 6));
		System.out.println(poly.evaluateAt(4, 6));
		System.out.println(poly.evaluateAt(5, 6));
		System.out.println(poly.evaluateAt(6, 6));
		System.out.println(poly.evaluateAt(7, 6));
		System.out.println(poly.evaluateAt(8, 6));

		System.out.println(poly.evaluateAt(1, 7));
		System.out.println(poly.evaluateAt(2, 7));
		System.out.println(poly.evaluateAt(3, 7));
		System.out.println(poly.evaluateAt(4, 7));
		System.out.println(poly.evaluateAt(5, 7));
		System.out.println(poly.evaluateAt(6, 7));
		System.out.println(poly.evaluateAt(7, 7));
		System.out.println(poly.evaluateAt(8, 7));

		System.out.println(poly.evaluateAt(1, 8));
		System.out.println(poly.evaluateAt(2, 8));
		System.out.println(poly.evaluateAt(3, 8));
		System.out.println(poly.evaluateAt(4, 8));
		System.out.println(poly.evaluateAt(5, 8));
		System.out.println(poly.evaluateAt(6, 8));
		System.out.println(poly.evaluateAt(7, 8));
		System.out.println(poly.evaluateAt(8, 8));

		System.out.println(poly.evaluateAt(1.5, 1.5));

		System.out.println(poly.evaluateAt(2.5, 2.5));

		System.out.println(poly.evaluateAt(3.5, 3.5));
		
		System.out.println(poly.evaluateAt(4.5, 4.5));

		System.out.println(poly.evaluateAt(5.5, 5.5));
		
		System.out.println(poly.evaluateAt(7.5, 7.5));
	}
}
