package com.wildbitsfoundry.etk4j.control;

public class Margins {
	public double GainMargin;
	public double PhaseMargin;
	public double GainCrossOverFrequency;
	public double PhaseCrossOverFrequency;
	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "Margins [GainMargin=" + GainMargin + ", PhaseMargin=" + PhaseMargin + ", GainCrossOverFrequency="
				+ GainCrossOverFrequency + ", PhaseCrossOverFrequency=" + PhaseCrossOverFrequency + "]";
	}
}
