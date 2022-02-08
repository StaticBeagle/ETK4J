package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.math.MathETK;

class Chebyshev2OrderCalculationStrategy implements FilterOrderCalculationStrategy {
    @Override
    public double calculateExactOrder(double nat, double gPass, double gStop) {
        return MathETK.acosh(Math.sqrt((gStop - 1.0) / (gPass - 1.0))) / MathETK.acosh(nat);
    }

    @Override
    public double calculateLowPassWn(int n, LowPassSpecs specs) {
        double wp = specs.getPassBandFrequency();
        double rp = specs.getPassBandRipple();
        double rs = specs.getStopBandAttenuation();

        double gStop = Math.pow(10, 0.1 * Math.abs(rs));
        double gPass = Math.pow(10, 0.1 * Math.abs(rp));
        double newFreq = Math.cosh(1.0 / n * MathETK.acosh(Math.sqrt((gStop - 1.0) / (gPass - 1.0))));
        double nat = wp * newFreq;
        return nat;
    }

    @Override
    public double calculateHighPassWn(int n, HighPassSpecs specs) {
        double wp = specs.getPassBandFrequency();
        double rp = specs.getPassBandRipple();
        double rs = specs.getStopBandAttenuation();

        double gStop = Math.pow(10, 0.1 * Math.abs(rs));
        double gPass = Math.pow(10, 0.1 * Math.abs(rp));
        double newFreq = Math.cosh(1.0 / n * MathETK.acosh(Math.sqrt((gStop - 1.0) / (gPass - 1.0))));
        double nat = wp / newFreq;
        return nat;
    }

    @Override
    public double[] calculateBandPassWn(int n, BandpassSpecs specs) {
        double wp1 = specs.getLowerPassBandFrequency();
        double wp2 = specs.getUpperPassBandFrequency();
        double rp = specs.getPassBandRipple();
        double rs = specs.getStopBandAttenuation();

        double gStop = Math.pow(10, 0.1 * Math.abs(rs));
        double gPass = Math.pow(10, 0.1 * Math.abs(rp));

        double newFreq = 1.0 / Math.cosh(1.0 / n * MathETK.acosh(Math.sqrt((gStop - 1.0) / (gPass - 1.0))));
        double[] nat = new double[2];
        nat[0] = (1.0 / (2.0 * newFreq) * (wp1 - wp2) +
                Math.sqrt(Math.pow(wp2 - wp1, 2) / (4.0 * Math.pow(newFreq, 2)) +
                wp2 * wp1));
        nat[1] = wp1 * wp2 / nat[0];
        return nat;
    }

    @Override
    public double[] calculateBandStopWn(int n, BandStopSpecs specs) {
        double wp1 = specs.getLowerPassBandFrequency();
        double wp2 = specs.getUpperPassBandFrequency();
        double rp = specs.getPassBandRipple();
        double rs = specs.getStopBandAttenuation();

        double gStop = Math.pow(10, 0.1 * Math.abs(rs));
        double gPass = Math.pow(10, 0.1 * Math.abs(rp));

        double newFreq = 1.0 / Math.cosh(1.0 / n * MathETK.acosh(Math.sqrt((gStop - 1.0) / (gPass - 1.0))));
        double[] nat = new double[2];
        nat[0] = newFreq / 2.0 * (wp1 - wp2) + Math.sqrt(Math.pow(newFreq, 2) * Math.pow((wp2 - wp1), 2) / 4.0 +
                wp2 * wp1);
        nat[1] = wp2 * wp1 / nat[0];
        return nat;
    }
}
