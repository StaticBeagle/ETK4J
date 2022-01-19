package com.wildbitsfoundry.etk4j.math.functions;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

/**
 * Interface describing a single variable complex valued function.
 */
public interface ComplexUnivariateFunction {
    Complex evaluateAt(Complex c);
}
