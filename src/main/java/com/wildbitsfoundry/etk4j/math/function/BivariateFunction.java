package com.wildbitsfoundry.etk4j.math.function;

/**
 * Interface describing a function of two variables.
 */
public interface BivariateFunction {
	double evaluateAt(double x, double y);
}
