package com.wildbitsfoundry.etk4j.control;

/**
 * The {@code Margins} class holds the calculation results
 * for the various performance metrics of a {@link TransferFunction}
 * i.e.
 * <ol>
 * <li>phase margin</li>
 * <li>gain margin</li>
 * <li>gain crossover frequency</li>
 * <li>phase crossover frequency</li>
 * </ol>
 */
public class Margins {
    private final double gainMargin;
    private final double phaseMargin;
    private final double phaseCrossoverFrequency;
    private final double gainCrossoverFrequency;

    public Margins(double gainMargin, double phaseMargin, double phaseCrossoverFrequency, double gainCrossoverFrequency) {
        this.gainMargin = gainMargin;
        this.phaseMargin = phaseMargin;
        this.phaseCrossoverFrequency = phaseCrossoverFrequency;
        this.gainCrossoverFrequency = gainCrossoverFrequency;
    }

    /**
     * Gain margin of the system.
     *
     * @return The gain margin of the system.
     */
    public double getGainMargin() {
        return gainMargin;
    }

    /**
     * Phase margin of the system.
     *
     * @return The phase margin of the system.
     */
    public double getPhaseMargin() {
        return phaseMargin;
    }

    /**
     * Phase crossover of the system. This is the frequency at which the gain margin is calculated.
     *
     * @return The phase crossover of the system.
     */
    public double getPhaseCrossoverFrequency() {
        return phaseCrossoverFrequency;
    }

    /**
     * Gain crossover of the system. This is the frequency at which the phase margin is calculated.
     *
     * @return The gain crossover for the system.
     */
    public double getGainCrossoverFrequency() {
        return gainCrossoverFrequency;
    }


    /* (non-Javadoc)
     * @see java.lang.Object#toString()
     */
    @Override
    public String toString() {
        return "Margins{" +
                "gainMargin=" + gainMargin +
                ", phaseMargin=" + phaseMargin +
                ", phaseCrossoverFrequency=" + phaseCrossoverFrequency +
                ", gainCrossoverFrequency=" + gainCrossoverFrequency +
                '}';
    }
}
