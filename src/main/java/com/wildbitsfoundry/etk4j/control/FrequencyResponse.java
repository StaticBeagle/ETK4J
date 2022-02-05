package com.wildbitsfoundry.etk4j.control;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

/**
 * The {@code Frequency Response} holds the {@link Complex }results of the evaluation of the
 * {@link TransferFunction} at a given set of frequencies.
 */
public class FrequencyResponse {
    private Complex[] response;
    private double[] w;

    FrequencyResponse(Complex[] response, double[] w) {
        this.response = response;
        this.w = w;
    }

    /**
     * Complex response of the system.
     * @return The {@link Complex} response of the system.
     */
    public Complex[] getResponse() {
        return response;
    }

    /**
     * Frequencies at which the system was evaluated.
     * @return An array of frequencies at which the system was evaluated.
     */
    public double[] getFrequencies() {
        return w;
    }
}