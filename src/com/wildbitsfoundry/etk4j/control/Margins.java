package com.wildbitsfoundry.etk4j.control;

public class Margins 
{
	private double[] _wcg; // Gain crossover frequency;
	private double[] _wcp; // Phase crossover frequency;
	private double _pm = Double.NaN;  // Phase margin;
	private double _gm = Double.NaN;  // Gain Margin;
	public Margins(final double[] wcg, final double[] wcp, double pm, double gm) 
	{		
		this._wcg = wcg;
		this._wcp = wcp;
		this._pm = pm;
		this._gm = gm;
	}

	public double getGainCrossOverFrequency()
	{
		return this._wcg[0];
	}
	
	public double[] getAllGainCrossOverFrequencies()
	{
		return this._wcg;
	}	
	
	public double getPhaseCrossOverFrequency()
	{
		return this._wcp[0];
	}
	
	
	public double[] getAllPhaseCrossOverFrequencies()
	{
		return this._wcp;
	}
	
	public double getPhaseMargin()
	{
		return this._pm;
	}
	
	public double getGainMargin()
	{
		return this._gm;
	}
	
	@Override
	public String toString()
	{
		StringBuilder sb = new StringBuilder();
		return sb.append("Gain Crossover [rad/sec]: ").append(String.format("%.4f", _wcg[0]))
				 .append(", Phase Crossover [rad/sec]: ").append(String.format("%.4f", _wcp[0]))
				 .append(", Phase Margin [deg]: ").append(String.format("%.4f", _pm))
				 .append(", Gain Margin [dB]: ").append(String.format("%.4f", _gm)).toString();
	}
}
