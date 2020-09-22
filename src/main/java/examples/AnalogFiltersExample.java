package examples;

import com.wildbitsfoundry.etk4j.signals.filters.AnalogFilter;
import com.wildbitsfoundry.etk4j.signals.filters.Butterworth;
import com.wildbitsfoundry.etk4j.signals.filters.Chebyshev;
import com.wildbitsfoundry.etk4j.signals.filters.Elliptic;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.LowPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.InverseChebyshev;

public class AnalogFiltersExample {

	public static void main(String[] args) {
		// Note on PassBandRipple (dB)
		// For Butterworth: Gain drop at the cutoff frequency
		// For Chebyshev & Elliptical: Ripple in the pass band
		// for Inverse Chebyshev: Ripple in the stop band

		// Specs for low pass filter
		LowPassSpecs lpSpecs = new LowPassSpecs();
		lpSpecs.setPassBandRipple(3.01); // 3.01 dB gain/ripple refer to note
		lpSpecs.setStopBandAttenuation(60.0); // 60 dB at the stop band
		lpSpecs.setPassBandFrequency(1.0); // 1 Hz cutoff frequency
		lpSpecs.setStopBandFrequency(10.0); // 10 Hz stop band frequency

		buildLowPassFilters(lpSpecs);

		// Specs for high pass filter
		HighPassSpecs hpSpecs = new HighPassSpecs();
		hpSpecs.setPassBandRipple(3.01); // 3.01 dB gain/ripple refer to note
		hpSpecs.setStopBandAttenuation(60.0); // 60 dB at the stop band
		hpSpecs.setPassBandFrequency(1.0); // 1 Hz cutoff frequency
		hpSpecs.setStopBandFrequency(0.1); // 0.1 Hz stop band frequency

		buildHighPassFilters(hpSpecs);

		// Specs for band pass filter
		BandPassSpecs bpSpecs = new BandPassSpecs();
		// The bandwidth of the filter starts at the LowerPassBandFrequency and
		// ends at the UpperPassBandFrequency. The filter has lower stop band
		// which is set LowerStopBandFrequency and the upper stop band can be set
		// with UpperStopBandFrequency. The attenuation at the stop bands can be
		// set with the LowerStopBandAttenuation and UpperStopBandAttenuation
		// respectively. In a frequency spectrum, the order of the frequencies will be:
		// LowerStopBandFrequency < LowerPassBandFrequency < UpperPassBandFrequency <
		// UpperStopBandFrequency
		bpSpecs.setLowerPassBandFrequency(190.0); // 190 Hz lower pass band frequency
		bpSpecs.setUpperPassBandFrequency(210.0); // 210 Hz upper pass band frequency
		bpSpecs.setLowerStopBandFrequency(180.0); // 180 Hz lower stop band frequency
		bpSpecs.setUpperStopBandFrequency(220.0); // 220 Hz upper stop band frequency
		bpSpecs.setPassBandRipple(0.2); // 0.2 dB gain/ripple refer to note
		bpSpecs.setLowerStopBandAttenuation(20.0); // 20 dB attenuation at the lower end of the skirt
		bpSpecs.setUpperStopBandAttenuation(20.0); // 20 dB attenuation at the upper end of the skirt

		buildBandPassFilters(bpSpecs);

		// Specs for low pass filter
		BandStopSpecs bsSpecs = new BandStopSpecs();
		// The notch of the filter starts at the LowerStopBandFrequency and
		// ends at the UpperStopBandFrequency. The filter has lower pass band
		// which is set LowerPassBandFrequency and the upper pass band can be set
		// with UpperPassBandFrequency. The attenuation at the notch can be
		// set with the StopBandAttenuation parameter and the attenuation/ripple
		// in the pass band can be set with the PassBandRipple parameter.
		// In a frequency spectrum, the order of the frequencies will be:
		// LowerPassBandFrequency < LowerStopBandFrequency < UpperStopBandFrequency <
		// UpperPassBandFrequency
		bsSpecs.setLowerPassBandFrequency(3.6e3); // 3600 Hz lower pass band frequency
		bsSpecs.setUpperPassBandFrequency(9.1e3); // 9100 Hz lower pass band frequency
		bsSpecs.setLowerStopBandFrequency(5.45e3); // 5450 Hz lower stop band frequency
		bsSpecs.setUpperStopBandFrequency(5.90e3); // 5900 Hz upper stop band frequency
		bsSpecs.setPassBandRipple(0.5); // 0.5 dB gain/ripple refer to note
		bsSpecs.setStopBandAttenuation(38.0); // 38 db attenuation at the notch

		buildBandStopFilters(bsSpecs);
		
		// Calculate the minimum order for a low pass Chebyshev to meet the following specs:
		// 0.2 dB ripple in the pass band
		// 60 dB attenuation in the stop band
		// 1000 Hz cutoff frequency
		// 10000 Hz stop band frequency
		
		// Step 1. Calculate the filter order required to meet the specs
		double ripple = 0.2;
		double attenuation = 60.0;
		double cutoff = 1000.0;
		double stopBand = 10000.0;
		int n = Chebyshev.getMinOrderNeeded(cutoff, stopBand, ripple, attenuation);
		System.out.println("Minimum order required to meet the specs is:");
		System.out.println(n);
	}

	public static void buildLowPassFilters(LowPassSpecs lpSpecs) {
		AnalogFilter bu = new Butterworth(lpSpecs);
		AnalogFilter cb1 = new Chebyshev(lpSpecs);
		AnalogFilter cb2 = new InverseChebyshev(lpSpecs);
		AnalogFilter el = new Elliptic(lpSpecs);
		
		System.out.println();
		System.out.println("//////////////////////////////////");
		System.out.println("//");
		System.out.println("// Low pass filter approximations");
		System.out.println("//");
		System.out.println("//////////////////////////////////");

		printTransferFunctions(bu, cb1, cb2, el);
	}

	public static void buildHighPassFilters(HighPassSpecs hpSpecs) {
		AnalogFilter bu = new Butterworth(hpSpecs);
		AnalogFilter cb1 = new Chebyshev(hpSpecs);
		AnalogFilter cb2 = new InverseChebyshev(hpSpecs);
		AnalogFilter el = new Elliptic(hpSpecs);
		
		System.out.println();
		System.out.println("//////////////////////////////////");
		System.out.println("//");
		System.out.println("// High pass filter approximations");
		System.out.println("//");
		System.out.println("//////////////////////////////////");

		printTransferFunctions(bu, cb1, cb2, el);
	}

	public static void buildBandPassFilters(BandPassSpecs bpSpecs) {
		AnalogFilter bu = new Butterworth(bpSpecs);
		AnalogFilter cb1 = new Chebyshev(bpSpecs);
		AnalogFilter cb2 = new InverseChebyshev(bpSpecs);
		AnalogFilter el = new Elliptic(bpSpecs);
		
		System.out.println();
		System.out.println("//////////////////////////////////");
		System.out.println("//");
		System.out.println("// Band pass filter approximations");
		System.out.println("//");
		System.out.println("//////////////////////////////////");

		printTransferFunctions(bu, cb1, cb2, el);
	}

	public static void buildBandStopFilters(BandStopSpecs bsSpecs) {
		AnalogFilter bu = new Butterworth(bsSpecs);
		AnalogFilter cb1 = new Chebyshev(bsSpecs);
		AnalogFilter cb2 = new InverseChebyshev(bsSpecs);
		AnalogFilter el = new Elliptic(bsSpecs);
		
		System.out.println();
		System.out.println("//////////////////////////////////");
		System.out.println("//");
		System.out.println("// Band stop filter approximations");
		System.out.println("//");
		System.out.println("//////////////////////////////////");

		printTransferFunctions(bu, cb1, cb2, el);
	}

	static void printTransferFunctions(AnalogFilter... analogFilters) {
		System.out.printf("Butterworth: %n%s%n%n", analogFilters[0]);
		System.out.printf("Chebyshev: %n%s%n%n", analogFilters[1]);
		System.out.printf("Inverse Chebyshev: %n%s%n%n", analogFilters[2]);
		System.out.printf("Elliptic: %n%s%n", analogFilters[3]);
	}
}
