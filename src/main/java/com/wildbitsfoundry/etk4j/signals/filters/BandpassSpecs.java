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
        if (lowerStopBandFrequency <= 0) {
            throw new IllegalArgumentException("The lower stop band frequency has to be greater than zero.");
        }
        if (lowerPassBandFrequency <= 0) {
            throw new IllegalArgumentException("The lower pass band frequency has to be greater than zero.");
        }
        if(upperPassBandFrequency <= 0) {
            throw new IllegalArgumentException("The upper pass band frequency has to be greater than zero.");
        }
        if(upperStopBandFrequency <= 0) {
            throw new IllegalArgumentException("The upper stop band frequency has to be greater than zero.");
        }
        if (lowerPassBandFrequency >= upperPassBandFrequency) {
            throw new IllegalArgumentException("The lower pass band frequency has to be less than the upper pass band frequency.");
        }
        if (lowerStopBandFrequency >= lowerPassBandFrequency) {
            throw new IllegalArgumentException("The lower stop band frequency has to be less than the lower pass band frequency.");
        }
        if (upperPassBandFrequency >= upperStopBandFrequency) {
            throw new IllegalArgumentException("The upper pass band frequency has to be less than the upper stop band frequency.");
        }
    }
}