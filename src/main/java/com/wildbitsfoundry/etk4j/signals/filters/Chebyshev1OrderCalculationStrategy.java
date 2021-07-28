package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.math.MathETK;

class Chebyshev1OrderCalculationStrategy implements FilterOrderCalculationStrategy {
    @Override
    public double calculateExactOrder(double nat, double gPass, double gStop) {
        return MathETK.acosh(Math.sqrt((gStop - 1.0) / (gPass - 1.0))) / MathETK.acosh(nat);
    }

    @Override
    public double calculateLowPassWn(int n, FilterSpecs.LowPassSpecs specs) {
        return specs.getPassBandFrequency();
    }

    @Override
    public double calculateHighPassWn(int n, FilterSpecs.HighPassSpecs specs) {
        return specs.getPassBandFrequency();
    }

    @Override
    public double[] calculateBandPassWn(int n, FilterSpecs.BandPassSpecs specs) {
        double[] wn = new double[2];
        wn[0] = specs.getLowerPassBandFrequency();
        wn[1] = specs.getUpperPassBandFrequency();
        return wn;
    }

    @Override
    public double[] calculateBandStopWn(int n, FilterSpecs.BandStopSpecs specs) {
        double wp1 = specs.getLowerPassBandFrequency();
        double wp2 = specs.getUpperPassBandFrequency();
        return new double[]{wp1, wp2};
    }
}
