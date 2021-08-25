package com.wildbitsfoundry.etk4j.math.functions;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

public interface ComplexUnivariateFunction {
    Complex evaluateAt(Complex c);
}
