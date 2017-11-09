package com.wildbitsfoundry.etk4j.systems.filters;

import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.BandPassSpecs;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.LowPassSpecs;

class ButterworthPrototype implements LowPassPrototype {

	@Override
	public AnalogFilter toHighPass(HighPassSpecs specs) {
		double fp = specs.PassBandFrequency;
		double fs = specs.StopBandFrequency;
		double ap = specs.PassBandAttenuation;
		double as = specs.StopBandAttenuation;
		
		return Butterworth.newHighPass(fp, fs, ap, as);
	}

	@Override
	public AnalogFilter toBandPass(BandPassSpecs specs) {
		double fp1 = specs.LowerPassBandFrequency;
		double fp2 = specs.UpperPassBandFrequency;
		double fs1 = specs.LowerStopBandFrequency;
		double fs2 = specs.UpperStopBandFrequency;
		double ap = specs.PassBandAttenuation;
		double as1 = specs.LowerStopBandAttenuation;
		double as2 = specs.UpperStopBandAttenuation;
		
		return Butterworth.newBandPass(fp1, fp2, fs1, fs2, ap, as1, as2);
	}

	@Override
	public AnalogFilter toBandStop(BandStopSpecs specs) {
		double fp1 = specs.LowerPassBandFrequency;
		double fp2 = specs.UpperPassBandFrequency;
		double fs1 = specs.LowerStopBandFrequency;
		double fs2 = specs.UpperStopBandFrequency;
		double amax = specs.PassBandAttenuation;
		double amin = specs.StopBandAttenuation;
		
		return Butterworth.newBandStop(fp1, fp2, fs1, fs2, amax, amin);
	}

	@Override
	public AnalogFilter buildLowPass(LowPassSpecs specs) {
		double fp = specs.PassBandFrequency;
		double fs = specs.StopBandFrequency;
		double ap = specs.PassBandAttenuation;
		double as = specs.StopBandAttenuation;
		return Butterworth.newLowPass(fp, fs, ap, as);
	}
}
