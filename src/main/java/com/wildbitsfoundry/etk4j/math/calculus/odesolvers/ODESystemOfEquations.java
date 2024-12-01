package com.wildbitsfoundry.etk4j.math.calculus.odesolvers;

public interface ODESystemOfEquations {
    double[] evaluateAt(double t, double[] y);
}
