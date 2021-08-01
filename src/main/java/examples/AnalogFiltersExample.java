package examples;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.signals.filters.*;
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
        lpSpecs.setPassBandRipple(1.5); // 1.5 dB gain/ripple refer to note
        lpSpecs.setStopBandAttenuation(60.0); // 60 dB at the stop band
        lpSpecs.setPassBandFrequency(2500); // 2500 Hz cutoff frequency
        lpSpecs.setStopBandFrequency(10000); // 10000 Hz stop band frequency

        buildLowPassFilters(lpSpecs);

        // Specs for high pass filter
        HighPassSpecs hpSpecs = new HighPassSpecs();
        hpSpecs.setPassBandRipple(0.2); // 0.2 dB gain/ripple refer to note
        hpSpecs.setStopBandAttenuation(60.0); // 60 dB at the stop band
        hpSpecs.setPassBandFrequency(12); // 12 Hz cutoff frequency
        hpSpecs.setStopBandFrequency(0.2); // 0.2 Hz stop band frequency

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
        bpSpecs.setStopBandAttenuation(20.0); // 20 dB attenuation in the stop band

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
        bsSpecs.setPassBandRipple(0.5); // 1.5 dB gain/ripple refer to note
        bsSpecs.setStopBandAttenuation(38.0); // 38 db attenuation at the notch

        buildBandStopFilters(bsSpecs);
    }

    public static void buildLowPassFilters(LowPassSpecs lpSpecs) {
        FilterOrderResults.OrderAndCutoffFrequency nW0 = ButterWorth.buttord(lpSpecs);
        TransferFunction bu = ButterWorth.newLowPass(nW0.getOrder(), nW0.getCutoffFrequency());

        nW0 = Chebyshev1.cheb1ord(lpSpecs);
        TransferFunction cb1 = Chebyshev1.newLowPass(nW0.getOrder(), lpSpecs.getPassBandRipple(), nW0.getCutoffFrequency());

        nW0 = Chebyshev2.cheb2ord(lpSpecs);
        TransferFunction cb2 = Chebyshev2.newLowPass(nW0.getOrder(), lpSpecs.getStopBandAttenuation(), nW0.getCutoffFrequency());

        nW0 = Elliptic.ellipord(lpSpecs);
        TransferFunction el = Elliptic.newLowPass(nW0.getOrder(), lpSpecs.getPassBandRipple(),
                lpSpecs.getStopBandAttenuation(), nW0.getCutoffFrequency());

        System.out.println();
        System.out.println("//////////////////////////////////");
        System.out.println("//");
        System.out.println("// Low pass filter approximations");
        System.out.println("//");
        System.out.println("//////////////////////////////////");

        printTransferFunctions(bu, cb1, cb2, el);
    }

    public static void buildHighPassFilters(HighPassSpecs hpSpecs) {
        FilterOrderResults.OrderAndCutoffFrequency nW0 = ButterWorth.buttord(hpSpecs);
        TransferFunction bu = ButterWorth.newHighPass(nW0.getOrder(), nW0.getCutoffFrequency());

        nW0 = Chebyshev1.cheb1ord(hpSpecs);
        TransferFunction cb1 = Chebyshev1.newHighPass(nW0.getOrder(), hpSpecs.getPassBandRipple(), nW0.getCutoffFrequency());

        nW0 = Chebyshev2.cheb2ord(hpSpecs);
        TransferFunction cb2 = Chebyshev2.newHighPass(nW0.getOrder(), hpSpecs.getStopBandAttenuation(), nW0.getCutoffFrequency());

        nW0 = Elliptic.ellipord(hpSpecs);
        TransferFunction el = Elliptic.newHighPass(nW0.getOrder(), hpSpecs.getPassBandRipple(),
                hpSpecs.getStopBandAttenuation(), nW0.getCutoffFrequency());

        System.out.println();
        System.out.println("//////////////////////////////////");
        System.out.println("//");
        System.out.println("// High pass filter approximations");
        System.out.println("//");
        System.out.println("//////////////////////////////////");

        printTransferFunctions(bu, cb1, cb2, el);
    }

    public static void buildBandPassFilters(BandPassSpecs bpSpecs) {
        FilterOrderResults.OrderAndCutoffFrequencies nW0W1 = ButterWorth.buttord(bpSpecs);
        TransferFunction bu = ButterWorth.newBandPass(nW0W1.getOrder(), nW0W1.getLowerCutoffFrequency(),
                nW0W1.getUpperCutoffFrequency());

        nW0W1 = Chebyshev1.cheb1ord(bpSpecs);
        TransferFunction cb1 = Chebyshev1.newBandPass(nW0W1.getOrder(), bpSpecs.getPassBandRipple(),
                nW0W1.getLowerCutoffFrequency(), nW0W1.getUpperCutoffFrequency());

        nW0W1 = Chebyshev2.cheb2ord(bpSpecs);
        TransferFunction cb2 = Chebyshev2.newBandPass(nW0W1.getOrder(), bpSpecs.getStopBandAttenuation(),
                nW0W1.getLowerCutoffFrequency(), nW0W1.getUpperCutoffFrequency());

        nW0W1 = Elliptic.ellipord(bpSpecs);
        TransferFunction el = Elliptic.newBandPass(nW0W1.getOrder(), bpSpecs.getPassBandRipple(),
                bpSpecs.getStopBandAttenuation(), nW0W1.getLowerCutoffFrequency(), nW0W1.getUpperCutoffFrequency());

        System.out.println();
        System.out.println("//////////////////////////////////");
        System.out.println("//");
        System.out.println("// Band pass filter approximations");
        System.out.println("//");
        System.out.println("//////////////////////////////////");

        printTransferFunctions(bu, cb1, cb2, el);
    }

    public static void buildBandStopFilters(BandStopSpecs bsSpecs) {
        FilterOrderResults.OrderAndCutoffFrequencies nW0W1 = ButterWorth.buttord(bsSpecs);
        TransferFunction bu = ButterWorth.newBandStop(nW0W1.getOrder(), nW0W1.getLowerCutoffFrequency(),
                nW0W1.getUpperCutoffFrequency());

        nW0W1 = Chebyshev1.cheb1ord(bsSpecs);
        TransferFunction cb1 = Chebyshev1.newBandStop(nW0W1.getOrder(), bsSpecs.getPassBandRipple(),
                nW0W1.getLowerCutoffFrequency(), nW0W1.getUpperCutoffFrequency());

        nW0W1 = Chebyshev2.cheb2ord(bsSpecs);
        TransferFunction cb2 = Chebyshev2.newBandStop(nW0W1.getOrder(), bsSpecs.getStopBandAttenuation(),
                nW0W1.getLowerCutoffFrequency(), nW0W1.getUpperCutoffFrequency());

        nW0W1 = Elliptic.ellipord(bsSpecs);
        TransferFunction el = Elliptic.newBandStop(nW0W1.getOrder(), bsSpecs.getPassBandRipple(),
                bsSpecs.getStopBandAttenuation(), nW0W1.getLowerCutoffFrequency(), nW0W1.getUpperCutoffFrequency());

        System.out.println();
        System.out.println("//////////////////////////////////");
        System.out.println("//");
        System.out.println("// Band stop filter approximations");
        System.out.println("//");
        System.out.println("//////////////////////////////////");

        printTransferFunctions(bu, cb1, cb2, el);
    }

    static void printTransferFunctions(TransferFunction... analogFilters) {
        System.out.printf("Butterworth: %n%s%n%n", analogFilters[0]);
        System.out.printf("Chebyshev: %n%s%n%n", analogFilters[1]);
        System.out.printf("Inverse Chebyshev: %n%s%n%n", analogFilters[2]);
        System.out.printf("Elliptic: %n%s%n", analogFilters[3]);
    }
}
