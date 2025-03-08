package com.wildbitsfoundry.etk4j.math.calculus.odesolver;

public interface OdeSystemOfEquations {
    double[] evaluateAt(double t, double[] y);
}
