package com.wildbitsfoundry.etk4j.signal.filters;

interface FilterOrderCalculationStrategy {
    double calculateExactOrder(double nat, double gPass, double gStop);
    double calculateLowPassWn(int n, LowPassSpecs specs);
    double calculateHighPassWn(int n, HighPassSpecs specs);
    double[] calculateBandPassWn(int n, BandpassSpecs specs);
    double[] calculateBandStopWn(int n, BandStopSpecs specs);

    default int calculateMinOrder(double nat, double gPass, double gStop) {
        return (int) Math.ceil(calculateExactOrder(nat, gPass, gStop));
    }
}
