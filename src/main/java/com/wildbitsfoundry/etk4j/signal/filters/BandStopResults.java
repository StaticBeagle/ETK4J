package com.wildbitsfoundry.etk4j.signal.filters;

public class BandStopResults {
    private int n;
    private double wn0;
    private double wn1;

    BandStopResults(int n, double wn0, double wn1) {
        this.n = n;
        this.wn0 = wn0;
        this.wn1 = wn1;
    }

    public int getOrder() {
        return n;
    }

    public double getLowerCutoffFrequency() {
        return wn0;
    }

    public double getUpperCutoffFrequency() {
        return wn1;
    }
}
