package com.wildbitsfoundry.etk4j.math.function;

/**
 * Interface describing a function of multiple variables.
 */
public interface MultivariateFunction {
	double evaluateAt(double... v);
}
