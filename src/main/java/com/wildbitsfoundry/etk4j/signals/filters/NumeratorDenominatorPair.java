package com.wildbitsfoundry.etk4j.signals.filters;

import java.util.Arrays;

public final class NumeratorDenominatorPair {

    private final double[] numerator;
    private final double[] denominator;

    NumeratorDenominatorPair(double[] numerator, double[] denominator) {
        this.numerator = numerator;
        this.denominator = denominator;
    }

    public double[] getNumerator() {
        return Arrays.copyOf(numerator, numerator.length);
    }

    public double[] getDenominator() {
        return Arrays.copyOf(denominator, denominator.length);
    }
}
