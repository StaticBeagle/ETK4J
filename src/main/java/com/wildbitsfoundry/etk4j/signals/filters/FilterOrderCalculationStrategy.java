package com.wildbitsfoundry.etk4j.signals.filters;

interface FilterOrderCalculationStrategy {
    double calculateExactOrder(double nat, double gPass, double gStop);
    double calculateLowPassWn(int n, FilterSpecs.LowPassSpecs specs);
    double calculateHighPassWn(int n, FilterSpecs.HighPassSpecs specs);
    double[] calculateBandPassWn(int n, FilterSpecs.BandPassSpecs specs);
    double[] calculateBandStopWn(int n, FilterSpecs.BandStopSpecs specs);

    default int calculateMinOrder(double nat, double gPass, double gStop) {
        return (int) Math.ceil(calculateExactOrder(nat, gPass, gStop));
    }
}
