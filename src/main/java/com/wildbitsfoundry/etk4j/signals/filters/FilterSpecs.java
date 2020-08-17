package com.wildbitsfoundry.etk4j.signals.filters;

public abstract class FilterSpecs {
	private static enum FilterType {
		LOWPASS,
		HIGHPASS,
		BANDPASS,
		BANDSTOP
	}

	protected FilterType _type;
	
	public static class LowPassSpecs extends FilterSpecs {
		public LowPassSpecs() {
			_type = FilterType.LOWPASS;
		}
		public double PassBandFrequency = 1;
		public double StopBandFrequency = 10;
		public double PassBandRipple = 0.2;
		public double StopBandAttenuation = 40;
	}
	
	public static class HighPassSpecs extends FilterSpecs {
		public HighPassSpecs() {
			_type = FilterType.HIGHPASS;
		}
		public double PassBandFrequency = 1;
		public double StopBandFrequency = 10;
		public double PassBandRipple = 0.2;
		public double StopBandAttenuation = 40;
	}
	
	public static class BandPassSpecs extends FilterSpecs {
		public BandPassSpecs() {
			_type = FilterType.BANDPASS;
		}
		public double LowerPassBandFrequency = 900;
		public double UpperPassBandFrequency = 1000;
		public double LowerStopBandFrequency = 90;
		public double UpperStopBandFrequency = 10000;
		public double PassBandRipple = 0.2;
		public double LowerStopBandAttenuation = 40;
		public double UpperStopBandAttenuation = 40;
	}
	
	public static class BandStopSpecs extends FilterSpecs {
		public BandStopSpecs() {
			_type = FilterType.BANDSTOP;
		}
		public double LowerPassBandFrequency = 3.6e3;
		public double UpperPassBandFrequency = 9.1e3;
		public double LowerStopBandFrequency = 5.45e3;
		public double UpperStopBandFrequency = 5.90e3;
		public double PassBandRipple = 1.5;
		public double StopBandAttenuation = 40;
	}

}
