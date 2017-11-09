package com.wildbitsfoundry.etk4j.systems.filters;

import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.BandPassSpecs;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.systems.filters.FilterSpecs.LowPassSpecs;

interface LowPassPrototype {
	AnalogFilter toHighPass(HighPassSpecs specs);
	AnalogFilter toBandPass(BandPassSpecs specs);
	AnalogFilter toBandStop(BandStopSpecs specs);
	AnalogFilter buildLowPass(LowPassSpecs specs);
}
