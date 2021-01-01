package com.wildbitsfoundry.etk4j.control;

/***
 * 
 * @author StaticBeagle
 * 
 * <p>
 * The {@code Margins} class holds the calculation results 
 * for the various performance metrics of a {@link TransferFunction} 
 * i.e.
 * <ol>
 * <li>phase margin</li>
 * <li>gain margin</li>
 * <li>gain crossover frequency</li>
 * <li>phase crossover frequency</li>
 * </ol>
 * </p>
 *
 */
public class Margins {
	private final double gainMargin;
	private final double phaseMargin;
	private final double gainCrossOverFrequency;
	private final double phaseCrossOverFrequency;

	public Margins(double gainMargin, double phaseMargin, double gainCrossOverFrequency, double phaseCrossOverFrequency) {
		this.gainMargin = gainMargin;
		this.phaseMargin = phaseMargin;
		this.gainCrossOverFrequency = gainCrossOverFrequency;
		this.phaseCrossOverFrequency = phaseCrossOverFrequency;
	}

	public double getGainMargin() {
		return gainMargin;
	}

	public double getPhaseMargin() {
		return phaseMargin;
	}

	public double getGainCrossOverFrequency() {
		return gainCrossOverFrequency;
	}

	public double getPhaseCrossOverFrequency() {
		return phaseCrossOverFrequency;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "Margins [GainMargin=" + gainMargin + ", PhaseMargin=" + phaseMargin + ", GainCrossOverFrequency="
				+ gainCrossOverFrequency + ", PhaseCrossOverFrequency=" + phaseCrossOverFrequency + "]";
	}
}
