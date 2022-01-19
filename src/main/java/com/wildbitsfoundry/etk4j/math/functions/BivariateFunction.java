package com.wildbitsfoundry.etk4j.math.functions;

/**
 * Interface describing a function of two variables.
 */
public interface BivariateFunction {
	double evaluateAt(double x, double y);
}
