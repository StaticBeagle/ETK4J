package com.wildbitsfoundry.etk4j.signal.filters;

/**
 * The {@code HighPassSpecs} class represents the design specifications for a high pass filter.
 */
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

    /**
     * Pass band (cutoff) frequency.
     * @param passBandFrequency The pass band frequency in radians per second (rad/s).
     */
    public void setPassBandFrequency(double passBandFrequency) {
        this.passBandFrequency = passBandFrequency;
    }

    public double getStopBandFrequency() {
        return stopBandFrequency;
    }

    /**
     * Stop band frequency.
     * @param stopBandFrequency The stop band frequency in radians per second (rad/s).
     */
    public void setStopBandFrequency(double stopBandFrequency) {
        this.stopBandFrequency = stopBandFrequency;
    }

    public double getPassBandRipple() {
        return passBandRipple;
    }

    /**
     * The pass band ripple.
     * @param passBandRipple Pass band ripple in decibels (dB). This value must be greater than zero.
     */
    public void setPassBandRipple(double passBandRipple) {
        this.passBandRipple = passBandRipple;
    }

    public double getStopBandAttenuation() {
        return stopBandAttenuation;
    }

    /**
     * The stop band attenuation.
     * @param stopBandAttenuation Stop band attenuation in decibels (dB). This value must be greater than zero.
     */
    public void setStopBandAttenuation(double stopBandAttenuation) {
        this.stopBandAttenuation = stopBandAttenuation;
    }

    void validate() {
        if(passBandFrequency <= 0) {
            throw new IllegalArgumentException("The pass band frequency must be grater than zero.");
        }
        if(stopBandFrequency <= 0) {
            throw new IllegalArgumentException("The stop band frequency must be grater than zero.");
        }
        if(passBandRipple <= 0) {
            throw new IllegalArgumentException("The pass band ripple must be grater than zero.");
        }
        if(stopBandAttenuation <= 0) {
            throw new IllegalArgumentException("The stop band attenuation must be grater than zero.");
        }
        if(stopBandAttenuation < passBandRipple) {
            throw new IllegalArgumentException("The stop band attenuation has to be greater than the pass band ripple");
        }
        if(stopBandFrequency >= passBandFrequency) {
            throw new IllegalArgumentException("The pass band frequency has to be greater than the stop band frequency");
        }
    }
}
