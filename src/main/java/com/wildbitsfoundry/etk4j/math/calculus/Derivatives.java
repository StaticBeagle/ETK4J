package com.wildbitsfoundry.etk4j.math.calculus;

import com.wildbitsfoundry.etk4j.math.functions.MultivariateFunction;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;

/**
 * The <code>Derivatives</code> class contains various methods that can be used for
 * computing the derivative of a function numerically.
 */
public final class Derivatives {

	private Derivatives() {
	}

	/**
	 * Forward difference.
	 * @param func Function to differentiate.
	 * @param x Argument at which to evaluate the derivative.
	 * @param h Step size of the differentiation.
	 * @return The derivative evaluated at <code>x</code> with step size <code>h</code>.
	 * @throws IllegalArgumentException If the step size is less than or equal to zero.
	 */
	public static double forwardDifference(UnivariateFunction func, double x, double h) {
		if(h <= 0) {
			throw new IllegalArgumentException("The step size h must be greater than zero");
		}
		return (func.evaluateAt(x + h) - func.evaluateAt(x)) / h;
	}

	/**
	 * Forward difference using a 3 point stencil.
	 * @param func Function to differentiate.
	 * @param x Argument at which to evaluate the derivative.
	 * @param h Step size of the differentiation.
	 * @return The derivative evaluated at <code>x</code> with step size <code>h</code>.
	 * @throws IllegalArgumentException If the step size is less than or equal to zero.
	 */
	public static double forwardDifference3Points(UnivariateFunction func, double x, double h) {
		if(h <= 0) {
			throw new IllegalArgumentException("The step size h must be greater than zero");
		}
		double num = -func.evaluateAt(x + 2 * h) + 4 * func.evaluateAt(x + h) - 3 * func.evaluateAt(x);
		return num / (2 * h);
	}

	/**
	 * Forward difference using a 5 point stencil.
	 * @param func Function to differentiate.
	 * @param x Argument at which to evaluate the derivative.
	 * @param h Step size of the differentiation.
	 * @return The derivative evaluated at <code>x</code> with step size <code>h</code>.
	 * @throws IllegalArgumentException If the step size is less than or equal to zero.
	 */
	public static double forwardDifference5Points(UnivariateFunction func, double x, double h) {
		if(h <= 0) {
			throw new IllegalArgumentException("The step size h must be greater than zero");
		}
		double num = -25 * func.evaluateAt(x) + 48 * func.evaluateAt(x + h) - 36 * func.evaluateAt(x + 2 * h)
				+ 16 * func.evaluateAt(x + 3 * h) - 3 * func.evaluateAt(x + 4 * h);
		return num / (12 * h);
	}

	/**
	 * Backward difference.
	 * @param func Function to differentiate.
	 * @param x Argument at which to evaluate the derivative.
	 * @param h Step size of the differentiation.
	 * @return The derivative evaluated at <code>x</code> with step size <code>h</code>.
	 * @throws IllegalArgumentException If the step size is less than or equal to zero.
	 */
	public static double backwardDifference(UnivariateFunction func, double x, double h) {
		if(h <= 0) {
			throw new IllegalArgumentException("The step size h must be greater than zero");
		}
		return (func.evaluateAt(x) - func.evaluateAt(x - h)) / h;
	}

	/**
	 * Backward difference using a 3 point stencil.
	 * @param func Function to differentiate.
	 * @param x Argument at which to evaluate the derivative.
	 * @param h Step size of the differentiation.
	 * @return The derivative evaluated at <code>x</code> with step size <code>h</code>.
	 * @throws IllegalArgumentException If the step size is less than or equal to zero.
	 */
	public static double backwardDifference3Points(UnivariateFunction func, double x, double h) {
		if(h <= 0) {
			throw new IllegalArgumentException("The step size h must be greater than zero");
		}
		double num = func.evaluateAt(x - 2 * h) - 4 * func.evaluateAt(x - h) + 3 * func.evaluateAt(x);
		return num / (2 * h);
	}

	/**
	 * Backward difference using a 5 point stencil.
	 * @param func Function to differentiate.
	 * @param x Argument at which to evaluate the derivative.
	 * @param h Step size of the differentiation.
	 * @return The derivative evaluated at <code>x</code> with step size <code>h</code>.
	 * @throws IllegalArgumentException if the step size is less than or equal to zero.
	 */
	public static double backwardDifference5Points(UnivariateFunction func, double x, double h) {
		if(h <= 0) {
			throw new IllegalArgumentException("The step size h must be greater than zero");
		}
		double num = 25 * func.evaluateAt(x) - 48 * func.evaluateAt(x - h) + 36 * func.evaluateAt(x - 2 * h)
				- 16 * func.evaluateAt(x - 3 * h) + 3 * func.evaluateAt(x - 4 * h);
		return num / (12 * h);
	}

	/**
	 * Centered difference.
	 * @param func Function to differentiate.
	 * @param x Argument at which to evaluate the derivative.
	 * @param h Step size of the differentiation.
	 * @return The derivative evaluated at <code>x</code> with step size <code>h</code>.
	 * @throws IllegalArgumentException If the step size is less than or equal to zero.
	 */
	public static double centeredDifference(UnivariateFunction func, double x, double h) {
		if(h <= 0) {
			throw new IllegalArgumentException("The step size h must be greater than zero");
		}
		return (func.evaluateAt(x + h) - func.evaluateAt(x - h)) / (2 * h);
	}

	/**
	 * Centered difference using a 5 point stencil.
	 * @param func Function to differentiate.
	 * @param x Argument at which to evaluate the derivative.
	 * @param h Step size of the differentiation.
	 * @return The derivative evaluated at <code>x</code> with step size <code>h</code>.
	 * @throws IllegalArgumentException If the step size is less than or equal to zero.
	 */
	public static double centeredDifference5Points(UnivariateFunction func, double x, double h) {
		if(h <= 0) {
			throw new IllegalArgumentException("The step size h must be greater than zero");
		}
		double num = 8 * (func.evaluateAt(x + h) - func.evaluateAt(x - h)) - func.evaluateAt(x + 2 * h)
				+ func.evaluateAt(x - 2 * h);
		return num / (12 * h);
	}

	/**
	 * Centered difference using a 7 point stencil.
	 * @param func Function to differentiate.
	 * @param x Argument at which to evaluate the derivative.
	 * @param h Step size of the differentiation.
	 * @return The derivative evaluated at <code>x</code> with step size <code>h</code>.
	 * @throws IllegalArgumentException If the step size is less than or equal to zero.
	 */
	public static double centeredDifference7Points(UnivariateFunction func, double x, double h) {
		if(h <= 0) {
			throw new IllegalArgumentException("The step size h must be greater than zero");
		}
		double num = 45 * (func.evaluateAt(x + h) - func.evaluateAt(x - h))
				+ 9 * (func.evaluateAt(x - 2 * h) - func.evaluateAt(x + 2 * h)) + func.evaluateAt(x + 3 * h)
				- func.evaluateAt(x - 3 * h);
		return num / (60 * h);
	}

	/**
	 * Ridder's algorithm for differentiation. <br>
	 * @see <a href="http://www.ibiblio.org/technicalc/tiplist/en/files/pdf/tips/tip6_26.pdf">Ridders</a>
	 * @param func Function to differentiate.
	 * @param x Argument at which to evaluate the derivative.
	 * @param h Step size of the differentiation.
	 * @return The derivative evaluated at <code>x</code> with step size <code>h</code>.
	 * @throws IllegalArgumentException If the step size is less than or equal to zero.
	 */
	public static double ridders(UnivariateFunction func, double x, double h) {
		if(h <= 0) {
			throw new IllegalArgumentException("The step size h must be greater than zero");
		}
		double con = 1.4;
		double con2 = con * con;
		int nTab = 10;
		double safe = 2.0;
		double error = Double.MAX_VALUE;
		double result = 0.0, errorTotal, factor;
		double[][] m = new double[nTab][nTab];

		m[0][0] = (func.evaluateAt(x + h) - func.evaluateAt(x - h)) / (2.0 * h);
		for(int i = 1; i < nTab; ++i) {
			h /= con;
			m[0][i] = (func.evaluateAt(x + h) - func.evaluateAt(x - h)) / (2.0 * h);
			factor = con2;
			for(int j = 1; j <= i; ++j) {
				m[j][i] = (m[j - 1][i] * factor - m[j - 1][i - 1]) / (factor - 1.0);
				factor *= con2;
				errorTotal = Math.max(Math.abs(m[j][i] - m[j - 1][i]), Math.abs(m[j][i] - m[j - 1][i - 1]));
				if(errorTotal <= error) {
					error = errorTotal;
					result = m[j][i];
				}
			}
			if(Math.abs(m[i][i] - m[i - 1][i - 1]) >= safe * error) {
				break;
			}
		}
		return result;
	}
}
