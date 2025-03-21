package com.wildbitsfoundry.etk4j.math.calculus;

import com.wildbitsfoundry.etk4j.math.function.MultivariateFunction;

public interface JacobianCalculationStrategy {
    double[][] calculateJacobian(MultivariateFunction[] functions, double[] x, Object... params);
}
