package com.wildbitsfoundry.etk4j.signals.filters;

public class BandpassSpecs {

    private double lowerPassBandFrequency;
    private double upperPassBandFrequency;
    private double lowerStopBandFrequency;
    private double upperStopBandFrequency;
    private double passBandRipple;
    private double stopBandAttenuation;

    public BandpassSpecs() {
    }

    public BandpassSpecs(BandpassSpecs specs) {
        this.lowerPassBandFrequency = specs.lowerPassBandFrequency;
        this.upperPassBandFrequency = specs.upperPassBandFrequency;
        this.lowerStopBandFrequency = specs.lowerStopBandFrequency;
        this.upperStopBandFrequency = specs.upperStopBandFrequency;
        this.passBandRipple = specs.passBandRipple;
        this.stopBandAttenuation = specs.stopBandAttenuation;
    }

    public double getLowerPassBandFrequency() {
        return lowerPassBandFrequency;
    }

    public void setLowerPassBandFrequency(double lowerPassBandFrequency) {
        this.lowerPassBandFrequency = lowerPassBandFrequency;
    }

    public double getUpperPassBandFrequency() {
        return upperPassBandFrequency;
    }

    public void setUpperPassBandFrequency(double upperPassBandFrequency) {
        this.upperPassBandFrequency = upperPassBandFrequency;
    }

    public double getLowerStopBandFrequency() {
        return lowerStopBandFrequency;
    }

    public void setLowerStopBandFrequency(double lowerStopBandFrequency) {
        this.lowerStopBandFrequency = lowerStopBandFrequency;
    }

    public double getUpperStopBandFrequency() {
        return upperStopBandFrequency;
    }

    public void setUpperStopBandFrequency(double upperStopBandFrequency) {
        this.upperStopBandFrequency = upperStopBandFrequency;
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