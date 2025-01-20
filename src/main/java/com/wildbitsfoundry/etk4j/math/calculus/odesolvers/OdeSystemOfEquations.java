package com.wildbitsfoundry.etk4j.math.calculus.odesolvers;

public interface OdeSystemOfEquations {
    double[] evaluateAt(double t, double[] y);
}
