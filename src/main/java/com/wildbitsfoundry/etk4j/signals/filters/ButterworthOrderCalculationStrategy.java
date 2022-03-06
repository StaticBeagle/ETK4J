package com.wildbitsfoundry.etk4j.signals.filters;

import java.util.Arrays;

class ButterworthOrderCalculationStrategy implements FilterOrderCalculationStrategy {
    @Override
    public double calculateExactOrder(double nat, double gPass, double gStop) {
        return Math.log10((gStop - 1.0) / (gPass - 1.0)) / (2 * Math.log10(nat));
    }

    @Override
    public double calculateLowPassWn(int n, LowPassSpecs specs) {
        double wp = specs.getPassBandFrequency();
        double w0 = calculateW0(n, specs.getPassBandRipple());;
        return w0 * wp;
    }

    @Override
    public double calculateHighPassWn(int n, HighPassSpecs specs) {
        double wp = specs.getPassBandFrequency();
        double w0 = calculateW0(n, specs.getPassBandRipple());;
        return wp / w0;
    }

    @Override
    public double[] calculateBandPassWn(int n, BandpassSpecs specs) {
        double wp1 = specs.getLowerPassBandFrequency();
        double wp2 = specs.getUpperPassBandFrequency();
        double w0 = calculateW0(n, specs.getPassBandRipple());
        double[] w0N = {-w0, w0};
        double[] wn = new double[2];

        wn[0] = (-w0N[0] * (wp2 - wp1) / 2.0 +
                Math.sqrt(w0N[0] * w0N[0] / 4.0 * Math.pow(wp2 - wp1, 2) +
                        wp1 * wp2));
        wn[1] = (-w0N[1] * (wp2 - wp1) / 2.0 +
                Math.sqrt(w0N[1] * w0N[1] / 4.0 * Math.pow(wp2 - wp1, 2) +
                        wp1 * wp2));
        wn[0] = Math.abs(wn[0]);
        wn[1] = Math.abs(wn[1]);

        Arrays.sort(wn);
        return wn;
    }

    @Override
    public double[] calculateBandStopWn(int n, BandStopSpecs specs) {
        double wp1 = specs.getLowerPassBandFrequency();
        double wp2 = specs.getUpperPassBandFrequency();
        double w0 = calculateW0(n, specs.getPassBandRipple());
        double[] wn = new double[2];
        double disc = Math.sqrt(Math.pow(wp2 - wp1, 2) + 4 * w0 * w0 * wp1 * wp2);
        wn[0] = ((wp2 - wp1) + disc) / (2 * w0);
        wn[1] = ((wp2 - wp1) - disc) / (2 * w0);
        wn[0] = Math.abs(wn[0]);
        wn[1] = Math.abs(wn[1]);
        return wn;
    }

    private static double calculateW0(int n, double rp) {
        if(n <= 0) {
            throw new IllegalArgumentException("The order of the Butterworth filter must be greater than zero.");
        }
        double gPass = Math.pow(10, 0.1 * Math.abs(rp));
        double w0 = Math.pow(gPass - 1.0, -1.0 / (2.0 * n));
        return w0;
    }
}
