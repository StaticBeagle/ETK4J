package com.wildbitsfoundry.etk4j.control;

/**
 * The {@code BodeResponse} class holds the magnitude and frequency response after evaluating the given
 * {@link TransferFunction} at an array of frequencies.
 */
public class BodeResponse {
    private double[] magnitudeIndB;
    private double[] phase;
    private double[] w;

    BodeResponse(double[] magnitudeIndB, double[] phase, double[] w) {
        this.magnitudeIndB = magnitudeIndB;
        TransferFunction.unwrapPhase(phase);
        this.phase = phase;
        this.w = w;
    }

    /**
     * Magnitude response of the system.
     * @return The magnitude response of the system. <br>
     * {@code Mag<sub>dB</sub> = 20 * log10(abs(response))}.
     */
    public double[] getMagnitudeIndB() {
        return magnitudeIndB;
    }

    /**
     * Phase response of the system.
     * @return The unwrapped phase response of the system is degrees.
     */
    public double[] getPhaseInDegrees() {
        return phase;
    }

    /**
     * Frequencies at which the system was evaluated.
     * @return An array of frequencies at which the system was evaluated.
     */
    public double[] getFrequencies() {
        return w;
    }
}