package com.wildbitsfoundry.etk4j.math.functions;

/**
 * Interface describing a function of multiple variables.
 */
public interface MultivariateFunction {
	double evaluateAt(double... v);
}
