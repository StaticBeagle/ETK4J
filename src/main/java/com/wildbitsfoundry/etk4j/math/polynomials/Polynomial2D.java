package com.wildbitsfoundry.etk4j.math.polynomials;

import java.util.Arrays;
import com.wildbitsfoundry.etk4j.math.functions.BivariateFunction;
import com.wildbitsfoundry.etk4j.math.linearalgebra.MatrixDense;
import com.wildbitsfoundry.etk4j.util.Grids;
import com.wildbitsfoundry.etk4j.util.DoubleArrays;

public class Polynomial2D implements BivariateFunction {
	private double[] _coefs;
	private int _n;
	private int _m;

	public Polynomial2D(double[] coefs, int n, int m) {
		// check that (n + 1) * (m + 1) = coefs.length
		_coefs = Arrays.copyOf(coefs, coefs.length);
		_n = n;
		_m = m;
	}

	/***
	 * Finds the coefficients of a polynomial P(x,y) of degree n + m that fits the
	 * data z best in a least-square sense. Always make sure that the number of
	 * points for x,y,z are the same and that the total number of points &ge; (n + 1)
	 * * (m + 1)
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
	 * @return Returns a bivariate polynomial of degree n + m fits the data z best
	 *         in a least-square sense
	 */
	public static Polynomial2D polyFit2D(double[] x, double[] y, double[] z, int n, int m) {
		double[][] vars = new double[2][];
		vars[0] = x;
		vars[1] = y;

		int[] nm = new int[] { n, m };
		double[] coefs = polyFitN(vars, z, nm);
		return new Polynomial2D(coefs, n, m);
	}

	/***
	 * Finds the coefficients of a polynomial P(x,y) of degree n + m that fits the
	 * data z best in a least-square sense. Always make sure that the number of
	 * points for x,y,z are the same and that the total number of points &ge; (n + 1)
	 * * (m + 1) in a regular grid
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
	 * @return Returns a bivariate polynomial of degree n + m fits the data z best
	 *         in a least-square sense
	 */
	public static Polynomial2D polyFit2D(double[] x, double[] y, double[][] z, int n, int m) {
		double[][] vars = new double[2][];
		// Build grid
		Grids.GridData gd = Grids.GridData.of(x, y);

		vars[0] = gd.X;
		vars[1] = gd.Y;

		int[] nm = new int[] { n, m };
		double[] coefficients = polyFitN(vars, DoubleArrays.flatten(z), nm);
		return new Polynomial2D(coefficients, n, m);
	}

	private static double[] polyFitN(double[][] vars, double[] z, int[] n) {

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
				vandermonde[i] = DoubleArrays.kronecker(vandermonde[i], tmp);
			}
		}

		MatrixDense A = new MatrixDense(vandermonde);
		MatrixDense b = new MatrixDense(z, z.length);

		MatrixDense c = A.solve(b);

		return c.getArray();
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
		StringBuilder result = new StringBuilder();
		for (int i = _m, k = 0; i >= 0; --i) {
			for (int j = _n; j >= 0; --j) {
				double coefAbs = Math.abs(_coefs[k]);
				if (coefAbs == 0.0) {
					++k;
					continue;
				}
				// Handle beginning
				String sign = "";
				if(k == 0) {
					if(_coefs[k] < 0) {
						sign = "-";
					} 
				} else {
					sign = " + ";
					if(_coefs[k] < 0) {
						sign = " - ";
					} 
				}
				// Handle end
				String coef = coefAbs == 1.0 ? sign : String.format("%s%.4g * ", sign, coefAbs);
				if (k == _coefs.length - 1) {
					coef = String.format("%s%.4g", sign, coefAbs);
				}
				result.append(coef);
				result.append(getVariablesAndPowers(i, j));
				++k;
			}
		}
		return result.toString();
	}
		
	private static String getVariablesAndPowers(int i, int j) {
		String x = null;
		String y = null;
		if (j == 0) {
			x = "";
		} else if (j == 1) {
			x = "x";
		} else {
			x = String.format("x^%d", j);
		}

		if (i == 0) {
			y = "";
		}
		else if (i == 1) {
			y = "y";
		} else {
			y = String.format("y^%d", i);
		}
		if (!x.isEmpty() && !y.isEmpty()) {
			return x + " * " + y;
		} else if (x.isEmpty()) {
			return y;
		} else {
			return x;
		}
	}
}
