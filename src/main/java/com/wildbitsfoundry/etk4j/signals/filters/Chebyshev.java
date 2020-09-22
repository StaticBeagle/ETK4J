package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.LowPassSpecs;

public class Chebyshev extends AnalogFilter {

	public Chebyshev(LowPassSpecs specs) {
		super(specs, new ChebyshevLowPassPrototype());
	}
	
	public Chebyshev(HighPassSpecs specs) {
		super(specs, new ChebyshevLowPassPrototype());
	}
	
	public Chebyshev(BandPassSpecs specs) {
		super(specs, new ChebyshevLowPassPrototype());
	}
	
	public Chebyshev(BandStopSpecs specs) {
		super(specs, new ChebyshevLowPassPrototype());
	}
	
	/***
	 * Calculate the minimum order required for Low-Pass Chebyshev filter
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
		return new ChebyshevLowPassPrototype().getMinOrderNeeded(fp, fs, ap, as);
	}
}
