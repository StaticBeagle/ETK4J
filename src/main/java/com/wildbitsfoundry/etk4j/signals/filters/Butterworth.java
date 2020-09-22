package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.LowPassSpecs;

public class Butterworth extends AnalogFilter {
	
	public Butterworth(LowPassSpecs specs) {
		super(specs, new ButterworthLowPassPrototype());
	}
	
	public Butterworth(HighPassSpecs specs) {
		super(specs, new ButterworthLowPassPrototype());
	}
	
	public Butterworth(BandPassSpecs specs) {
		super(specs, new ButterworthLowPassPrototype());
	}
	
	public Butterworth(BandStopSpecs specs) {
		super(specs, new ButterworthLowPassPrototype());
	}
	

	/***
	 * Calculate the minimum order required for Low-Pass Butterworth filter
	 * 
	 * @param fp
	 *            passband frequency in Hertz
	 * @param fs
	 *            stopband frequency in Hertz
	 * @param ap
	 *            passband attenuation in dB
	 * @param as
	 *            stopband attenuation in dB
	 * @return
	 */
	public static int getMinOrderNeeded(double fp, double fs, double ap, double as) {
		return new ButterworthLowPassPrototype().getMinOrderNeeded(fp, fs, ap, as);
	}
}
