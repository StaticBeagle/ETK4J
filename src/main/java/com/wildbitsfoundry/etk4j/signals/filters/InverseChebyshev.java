package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.LowPassSpecs;

public class InverseChebyshev extends AnalogFilter {

	public InverseChebyshev(LowPassSpecs specs) {
		super(specs, new InverseChebyshevLowPassPrototype());
	}
	
	public InverseChebyshev(HighPassSpecs specs) {
		super(specs, new InverseChebyshevLowPassPrototype());
	}
	
	public InverseChebyshev(BandPassSpecs specs) {
		super(specs, new InverseChebyshevLowPassPrototype());
	}
	
	public InverseChebyshev(BandStopSpecs specs) {
		super(specs, new InverseChebyshevLowPassPrototype());
	}
	
	/***
	 * Calculate the minimum order required for Low-Pass Inverse Chebyshev filter
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
		return new InverseChebyshevLowPassPrototype().getMinOrderNeeded(fp, fs, ap, as);
	}
}
