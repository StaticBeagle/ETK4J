package com.wildbitsfoundry.etk4j.signals.filters;

public class LowPassSpecs {

    private double passBandFrequency;
    private double stopBandFrequency;
    private double passBandRipple;
    private double stopBandAttenuation;

    public LowPassSpecs() {}

    public LowPassSpecs(LowPassSpecs specs) {
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

    void validate() {
        if (passBandRipple < 0) {
            throw new IllegalArgumentException("The pass band ripple cannot be less than zero.");
        }
        if (stopBandAttenuation < 0) {
            throw new IllegalArgumentException("The stop band attenuation cannot be less than zero.");
        }
        if(passBandRipple >= stopBandAttenuation) {
            throw new IllegalArgumentException("The stop band attenuation has to be greater than the pass band ripple.");
        }
    }
}