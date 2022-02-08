package com.wildbitsfoundry.etk4j.signals.filters;

public class HighPassSpecs {

    private double passBandFrequency;
    private double stopBandFrequency;
    private double passBandRipple;
    private double stopBandAttenuation;

    public HighPassSpecs() {}

    public HighPassSpecs(HighPassSpecs specs) {
        this.passBandFrequency = specs.passBandFrequency;
        this.stopBandFrequency = specs.stopBandFrequency;
        this.passBandRipple = specs.passBandRipple;
        this.stopBandAttenuation = specs.stopBandAttenuation;
    }

    public double getPassBandFrequency() {
        return passBandFrequency;
    }

    public void setPassBandFrequency(double passBandFrequency) {
        this.passBandFrequency = passBandFrequency;
    }

    public double getStopBandFrequency() {
        return stopBandFrequency;
    }

    public void setStopBandFrequency(double stopBandFrequency) {
        this.stopBandFrequency = stopBandFrequency;
    }

    public double getPassBandRipple() {
        return passBandRipple;
    }

    public void setPassBandRipple(double passBandRipple) {
        this.passBandRipple = passBandRipple;
    }

    public double getStopBandAttenuation() {
        return stopBandAttenuation;
    }

    public void setStopBandAttenuation(double stopBandAttenuation) {
        this.stopBandAttenuation = stopBandAttenuation;
    }
}
