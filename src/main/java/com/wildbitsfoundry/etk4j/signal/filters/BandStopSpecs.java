package com.wildbitsfoundry.etk4j.signal.filters;

/**
 * The {@code BandStopSpec} class represents the design specifications for a band stop (notch) filter.
 */
public class BandStopSpecs {

    private double lowerPassBandFrequency;
    private double upperPassBandFrequency;
    private double lowerStopBandFrequency;
    private double upperStopBandFrequency;
    private double passBandRipple;
    private double stopBandAttenuation;

    public BandStopSpecs() {}

    public BandStopSpecs(BandStopSpecs specs) {
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

    /**
     * Lower pass band frequency.
     * @param lowerPassBandFrequency The lower pass band frequency in radians per second (rad/s).
     */
    public void setLowerPassBandFrequency(double lowerPassBandFrequency) {
        this.lowerPassBandFrequency = lowerPassBandFrequency;
    }

    public double getUpperPassBandFrequency() {
        return upperPassBandFrequency;
    }

    /**
     * Upper pass band frequency.
     * @param upperPassBandFrequency The upper pass band frequency in radians per second (rad/s).
     */
    public void setUpperPassBandFrequency(double upperPassBandFrequency) {
        this.upperPassBandFrequency = upperPassBandFrequency;
    }

    public double getLowerStopBandFrequency() {
        return lowerStopBandFrequency;
    }

    /**
     * Lower stop band frequency.
     * @param lowerStopBandFrequency The lower stop band frequency in radians per second (rad/s).
     */
    public void setLowerStopBandFrequency(double lowerStopBandFrequency) {
        this.lowerStopBandFrequency = lowerStopBandFrequency;
    }

    public double getUpperStopBandFrequency() {
        return upperStopBandFrequency;
    }

    /**
     * Upper stop band frequency.
     * @param upperStopBandFrequency The upper stop band frequency in radians per second (rad/s).
     */
    public void setUpperStopBandFrequency(double upperStopBandFrequency) {
        this.upperStopBandFrequency = upperStopBandFrequency;
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
        if (passBandRipple <= 0) {
            throw new IllegalArgumentException("The pass band ripple has to be greater than zero.");
        }
        if (stopBandAttenuation <= 0) {
            throw new IllegalArgumentException("The stop band attenuation has to be greater than zero.");
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
        if (lowerStopBandFrequency >= upperStopBandFrequency) {
            throw new IllegalArgumentException("The lower stop band frequency has to be less than the upper stop band frequency.");
        }
        if (lowerPassBandFrequency >= lowerStopBandFrequency) {
            throw new IllegalArgumentException("The lower pass band frequency has to be less than the lower stop band frequency.");
        }
        if (upperStopBandFrequency >= upperPassBandFrequency) {
            throw new IllegalArgumentException("The upper stop band frequency has to be less than the upper pass band frequency.");
        }
    }
}