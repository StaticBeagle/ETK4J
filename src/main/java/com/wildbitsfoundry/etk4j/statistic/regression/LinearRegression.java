package com.wildbitsfoundry.etk4j.statistic.regression;

/**
 * The {@code LinearRegression} class implements a liner fit in the least square sense for a given set of (x, y) points.
 */
public class LinearRegression extends PolynomialRegression {

	public LinearRegression(double[] x, double[] y) {
		super(x, y, 1);
	}
}
