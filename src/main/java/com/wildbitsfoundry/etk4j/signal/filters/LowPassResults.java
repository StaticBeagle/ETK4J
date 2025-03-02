package com.wildbitsfoundry.etk4j.signal.filters;

public class LowPassResults {
    private int n;
    private double wn;

    LowPassResults(int n, double wn) {
        this.n = n;
        this.wn = wn;
    }

    public int getOrder() {
        return n;
    }

    public double getCutoffFrequency() {
        return wn;
    }
}
