package com.wildbitsfoundry.etk4j.systems.continuoustime.filters;

import com.wildbitsfoundry.etk4j.systems.continuoustime.filters.FilterSpecs.BandPassSpecs;
import com.wildbitsfoundry.etk4j.systems.continuoustime.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.systems.continuoustime.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.systems.continuoustime.filters.FilterSpecs.LowPassSpecs;

interface LowPassPrototype {
	AnalogFilter toHighPass(HighPassSpecs specs);
	AnalogFilter toBandPass(BandPassSpecs specs);
	AnalogFilter toBandStop(BandStopSpecs specs);
	AnalogFilter buildLowPass(LowPassSpecs specs);
}
