package examples;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.signals.filters.AnalogFilter;
import com.wildbitsfoundry.etk4j.signals.filters.ApproximationType;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.LowPassSpecs;

public class AnalogFiltersExample {

	public static void main(String[] args) {
		// Note on PassBandRipple (dB)
		// For Butterworth: Gain drop at the cutoff frequency
		// For Chebyshev & Elliptical: Ripple in the pass band
		// for Inverse Chebyshev: Ripple in the stop band

		// Specs for low pass filter
		LowPassSpecs lpSpecs = new LowPassSpecs();
		lpSpecs.PassBandRipple = 3.01; // 3.01 dB gain/ripple refer to note
		lpSpecs.StopBandAttenuation = 60; // 60 dB at the stop band
		lpSpecs.PassBandFrequency = 1.0; // 1 Hz cutoff frequency
		lpSpecs.StopBandFrequency = 10.0; // 10 Hz stop band frequency

		buildLowPassFilters(lpSpecs);

		// Specs for high pass filter
		HighPassSpecs hpSpecs = new HighPassSpecs();
		hpSpecs.PassBandRipple = 3.01; // 3.01 dB gain/ripple refer to note
		hpSpecs.StopBandAttenuation = 60; // 60 dB at the stop band
		hpSpecs.PassBandFrequency = 1.0; // 1 Hz cutoff frequency
		hpSpecs.StopBandFrequency = 0.1; // 0.1 Hz stop band frequency

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
		bpSpecs.LowerPassBandFrequency = 190; // 190 Hz lower pass band frequency
		bpSpecs.UpperPassBandFrequency = 210; // 210 Hz upper pass band frequency
		bpSpecs.LowerStopBandFrequency = 180; // 180 Hz lower stop band frequency
		bpSpecs.UpperStopBandFrequency = 220; // 220 Hz upper stop band frequency
		bpSpecs.PassBandRipple = 0.2; // 0.2 dB gain/ripple refer to note
		bpSpecs.LowerStopBandAttenuation = 20; // 20 dB attenuation at the lower end of the skirt
		bpSpecs.UpperStopBandAttenuation = 20; // 20 dB attenuation at the upper end of the skirt

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
		bsSpecs.LowerPassBandFrequency = 3.6e3; // 3600 Hz lower pass band frequency
		bsSpecs.UpperPassBandFrequency = 9.1e3; // 9100 Hz lower pass band frequency
		bsSpecs.LowerStopBandFrequency = 5.45e3; // 5450 Hz lower stop band frequency
		bsSpecs.UpperStopBandFrequency = 5.90e3; // 5900 Hz upper stop band frequency
		bsSpecs.PassBandRipple = 0.5; // 0.5 dB gain/ripple refer to note
		bsSpecs.StopBandAttenuation = 38; // 38 db attenuation at the notch

		buildBandStopFilters(bsSpecs);
		
		// Design a low pass Chebyshev to meet the following specs:
		// 0.2 dB ripple in the pass band
		// 60 dB attenuation in the stop band
		// 1000 Hz cutoff frequency
		// 10000 Hz stop band frequency
		
		// Step 1. Calculate the filter order required to meet the specs
		double ripple = 0.2;
		double attenuation = 60;
		double cutoff = 1000;
		double stopBand = 10000;
		ApproximationType type = ApproximationType.CHEBYSHEV;
		int n = AnalogFilter.getMinOrderNeeded(cutoff, stopBand, ripple, attenuation, type);
		// Step 2. Calculate the filter approximation
		AnalogFilter filter = AnalogFilter.newLowPass(n, ripple, attenuation, type);
		
		// Step 3. Print the filter coefficients
		System.out.println();
		System.out.println("Filter numerator:");
		System.out.println(Arrays.toString(filter.getNumerator()));
		System.out.println("Filter denominator:");
		System.out.println(Arrays.toString(filter.getDenominator()));
		System.out.println("Filter transfer function:");
		System.out.println(filter);
	}

	public static void buildLowPassFilters(LowPassSpecs lpSpecs) {
		// Build Butterworth approximation
		AnalogFilter bu = AnalogFilter.newLowPass(lpSpecs, ApproximationType.BUTTERWORTH);
		AnalogFilter cb1 = AnalogFilter.newLowPass(lpSpecs, ApproximationType.CHEBYSHEV);
		AnalogFilter cb2 = AnalogFilter.newLowPass(lpSpecs, ApproximationType.INVERSE_CHEBYSHEV);
		AnalogFilter el = AnalogFilter.newLowPass(lpSpecs, ApproximationType.ELLIPTIC);

		printTransferFunctions(bu, cb1, cb2, el);
	}

	public static void buildHighPassFilters(HighPassSpecs hpSpecs) {
		// Build Butterworth approximation
		AnalogFilter bu = AnalogFilter.newHighPass(hpSpecs, ApproximationType.BUTTERWORTH);
		AnalogFilter cb1 = AnalogFilter.newHighPass(hpSpecs, ApproximationType.CHEBYSHEV);
		AnalogFilter cb2 = AnalogFilter.newHighPass(hpSpecs, ApproximationType.INVERSE_CHEBYSHEV);
		AnalogFilter el = AnalogFilter.newHighPass(hpSpecs, ApproximationType.ELLIPTIC);

		printTransferFunctions(bu, cb1, cb2, el);
	}

	public static void buildBandPassFilters(BandPassSpecs bpSpecs) {
		// Build Butterworth approximation
		AnalogFilter bu = AnalogFilter.newBandPass(bpSpecs, ApproximationType.BUTTERWORTH);
		AnalogFilter cb1 = AnalogFilter.newBandPass(bpSpecs, ApproximationType.CHEBYSHEV);
		AnalogFilter cb2 = AnalogFilter.newBandPass(bpSpecs, ApproximationType.INVERSE_CHEBYSHEV);
		AnalogFilter el = AnalogFilter.newBandPass(bpSpecs, ApproximationType.ELLIPTIC);

		printTransferFunctions(bu, cb1, cb2, el);
	}

	public static void buildBandStopFilters(BandStopSpecs bsSpecs) {
		// Build Butterworth approximation
		AnalogFilter bu = AnalogFilter.newBandStop(bsSpecs, ApproximationType.BUTTERWORTH);
		AnalogFilter cb1 = AnalogFilter.newBandStop(bsSpecs, ApproximationType.CHEBYSHEV);
		AnalogFilter cb2 = AnalogFilter.newBandStop(bsSpecs, ApproximationType.INVERSE_CHEBYSHEV);
		AnalogFilter el = AnalogFilter.newBandStop(bsSpecs, ApproximationType.ELLIPTIC);

		printTransferFunctions(bu, cb1, cb2, el);
	}

	static void printTransferFunctions(AnalogFilter... analogFilters) {
		System.out.println();
		System.out.println("Band stop filter approximations");
		System.out.println("------------------------------");
		System.out.printf("Butterworth: %n%s%n%n", analogFilters[0]);
		System.out.printf("Chebyshev: %n%s%n%n", analogFilters[1]);
		System.out.printf("Inverse Chebyshev: %n%s%n%n", analogFilters[2]);
		System.out.printf("Elliptic: %n%s%n", analogFilters[3]);
	}
}
