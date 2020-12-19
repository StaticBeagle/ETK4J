package examples;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
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

        System.out.println();
        double[] b = {1.0, 2.0, 3.0, 4.0, 5.0};
        double[] a = {5.0, 4.0, 3.0, 2.0, 1.0};
        TransferFunction tf = AnalogFilter.lpTobs(b, a, 2, 3);
        tf.normalize();
        System.out.println(tf);
    }

    public static void buildLowPassFilters(LowPassSpecs lpSpecs) {
        AnalogFilter bu = AnalogFilter.newLowPass(lpSpecs, ApproximationType.BUTTERWORTH);
        AnalogFilter cb1 = AnalogFilter.newLowPass(lpSpecs, ApproximationType.CHEBYSHEV);
        AnalogFilter cb2 = AnalogFilter.newLowPass(lpSpecs, ApproximationType.INVERSE_CHEBYSHEV);
        AnalogFilter el = AnalogFilter.newLowPass(lpSpecs, ApproximationType.ELLIPTIC);

        System.out.println();
        System.out.println("//////////////////////////////////");
        System.out.println("//");
        System.out.println("// Low pass filter approximations");
        System.out.println("//");
        System.out.println("//////////////////////////////////");

        printTransferFunctions(bu, cb1, cb2, el);
    }

    public static void buildHighPassFilters(HighPassSpecs hpSpecs) {
        AnalogFilter bu = AnalogFilter.newHighPass(hpSpecs, ApproximationType.BUTTERWORTH);
        AnalogFilter cb1 = AnalogFilter.newHighPass(hpSpecs, ApproximationType.CHEBYSHEV);
        AnalogFilter cb2 = AnalogFilter.newHighPass(hpSpecs, ApproximationType.INVERSE_CHEBYSHEV);
        AnalogFilter el = AnalogFilter.newHighPass(hpSpecs, ApproximationType.ELLIPTIC);

        System.out.println();
        System.out.println("//////////////////////////////////");
        System.out.println("//");
        System.out.println("// High pass filter approximations");
        System.out.println("//");
        System.out.println("//////////////////////////////////");

        printTransferFunctions(bu, cb1, cb2, el);
    }

    public static void buildBandPassFilters(BandPassSpecs bpSpecs) {
        AnalogFilter bu = AnalogFilter.newBandPass(bpSpecs, ApproximationType.BUTTERWORTH);
        AnalogFilter cb1 = AnalogFilter.newBandPass(bpSpecs, ApproximationType.CHEBYSHEV);
        AnalogFilter cb2 = AnalogFilter.newBandPass(bpSpecs, ApproximationType.INVERSE_CHEBYSHEV);
        AnalogFilter el = AnalogFilter.newBandPass(bpSpecs, ApproximationType.ELLIPTIC);


        System.out.println();
        System.out.println("//////////////////////////////////");
        System.out.println("//");
        System.out.println("// Band pass filter approximations");
        System.out.println("//");
        System.out.println("//////////////////////////////////");

        printTransferFunctions(bu, cb1, cb2, el);
    }

    public static void buildBandStopFilters(BandStopSpecs bsSpecs) {
        AnalogFilter bu = AnalogFilter.newBandStop(bsSpecs, ApproximationType.BUTTERWORTH);
        System.out.println(Arrays.toString(bu.getNumerator()));
        System.out.println(Arrays.toString(bu.getDenominator()));
        AnalogFilter cb1 = AnalogFilter.newBandStop(bsSpecs, ApproximationType.CHEBYSHEV);
        System.out.println(Arrays.toString(cb1.getNumerator()));
        System.out.println(Arrays.toString(cb1.getDenominator()));
        AnalogFilter cb2 = AnalogFilter.newBandStop(bsSpecs, ApproximationType.INVERSE_CHEBYSHEV);
        System.out.println(Arrays.toString(cb2.getNumerator()));
        System.out.println(Arrays.toString(cb2.getDenominator()));
        AnalogFilter el = AnalogFilter.newBandStop(bsSpecs, ApproximationType.ELLIPTIC);
        System.out.println(Arrays.toString(el.getNumerator()));
        System.out.println(Arrays.toString(el.getDenominator()));

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
