package com.wildbitsfoundry.etk4j.signals.filters;

import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.LowPassSpecs;

public class Elliptic extends AnalogFilter {

	public Elliptic(LowPassSpecs specs) {
		super(specs, new EllipticLowPassPrototype());
	}
	
	public Elliptic(HighPassSpecs specs) {
		super(specs, new EllipticLowPassPrototype());
	}
	
	public Elliptic(BandPassSpecs specs) {
		super(specs, new EllipticLowPassPrototype());
	}
	
	public Elliptic(BandStopSpecs specs) {
		super(specs, new EllipticLowPassPrototype());
	}
	
	/***
	 * Calculate the minimum order required for Low-Pass Elliptic filter
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
		return new EllipticLowPassPrototype().getMinOrderNeeded(fp, fs, ap, as);
	}
}
