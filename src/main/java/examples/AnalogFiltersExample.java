package examples;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.signals.filters.*;

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
        lpSpecs.setPassBandFrequency(2500); // 2500 rad/s cutoff frequency
        lpSpecs.setStopBandFrequency(10000); // 10000 rad/s stop band frequency

        buildLowPassFilters(lpSpecs);

        // Specs for high pass filter
        HighPassSpecs hpSpecs = new HighPassSpecs();
        hpSpecs.setPassBandRipple(0.2); // 0.2 dB gain/ripple refer to note
        hpSpecs.setStopBandAttenuation(60.0); // 60 dB at the stop band
        hpSpecs.setPassBandFrequency(12); // 12 rad/s cutoff frequency
        hpSpecs.setStopBandFrequency(0.2); // 0.2 rad/s stop band frequency

        buildHighPassFilters(hpSpecs);

        // Specs for band pass filter
        BandpassSpecs bpSpecs = new BandpassSpecs();
        // The bandwidth of the filter starts at the LowerPassBandFrequency and
        // ends at the UpperPassBandFrequency. The filter has lower stop band
        // which is set LowerStopBandFrequency and the upper stop band can be set
        // with UpperStopBandFrequency.
        // In a frequency spectrum, the order of the frequencies will be:
        // LowerStopBandFrequency < LowerPassBandFrequency < UpperPassBandFrequency < UpperStopBandFrequency
        bpSpecs.setLowerPassBandFrequency(190.0); // 190 rad/s lower pass band frequency
        bpSpecs.setUpperPassBandFrequency(210.0); // 210 rad/s upper pass band frequency
        bpSpecs.setLowerStopBandFrequency(180.0); // 180 rad/s lower stop band frequency
        bpSpecs.setUpperStopBandFrequency(220.0); // 220 rad/s upper stop band frequency
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
        // LowerPassBandFrequency < LowerStopBandFrequency < UpperStopBandFrequency < UpperPassBandFrequency
        bsSpecs.setLowerPassBandFrequency(3.6e3); // 3600 rad/s lower pass band frequency
        bsSpecs.setUpperPassBandFrequency(9.1e3); // 9100 rad/s lower pass band frequency
        bsSpecs.setLowerStopBandFrequency(5.45e3); // 5450 rad/s lower stop band frequency
        bsSpecs.setUpperStopBandFrequency(5.90e3); // 5900 rad/s upper stop band frequency
        bsSpecs.setPassBandRipple(0.5); // 1.5 dB gain/ripple refer to note
        bsSpecs.setStopBandAttenuation(38.0); // 38 db attenuation at the notch

        buildBandStopFilters(bsSpecs);
    }

    public static void buildLowPassFilters(LowPassSpecs lpSpecs) {
        LowPassResults lpr = ButterWorth.buttOrd(lpSpecs);
        TransferFunction bu = ButterWorth.newLowPass(lpr.getOrder(), lpr.getCutoffFrequency());

        lpr = Chebyshev1.cheb1Ord(lpSpecs);
        TransferFunction cb1 = Chebyshev1.newLowPass(lpr.getOrder(), lpSpecs.getPassBandRipple(), lpr.getCutoffFrequency());

        lpr = Chebyshev2.cheb2Ord(lpSpecs);
        TransferFunction cb2 = Chebyshev2.newLowPass(lpr.getOrder(), lpSpecs.getStopBandAttenuation(), lpr.getCutoffFrequency());

        lpr = Elliptic.ellipOrd(lpSpecs);
        TransferFunction el = Elliptic.newLowPass(lpr.getOrder(), lpSpecs.getPassBandRipple(),
                lpSpecs.getStopBandAttenuation(), lpr.getCutoffFrequency());

        System.out.println();
        System.out.println("//////////////////////////////////");
        System.out.println("//");
        System.out.println("// Low pass filter approximations");
        System.out.println("//");
        System.out.println("//////////////////////////////////");

        printTransferFunctions(bu, cb1, cb2, el);
    }

    public static void buildHighPassFilters(HighPassSpecs hpSpecs) {
        HighPassResults hpr = ButterWorth.buttOrd(hpSpecs);
        TransferFunction bu = ButterWorth.newHighPass(hpr.getOrder(), hpr.getCutoffFrequency());

        hpr = Chebyshev1.cheb1Ord(hpSpecs);
        TransferFunction cb1 = Chebyshev1.newHighPass(hpr.getOrder(), hpSpecs.getPassBandRipple(), hpr.getCutoffFrequency());

        hpr = Chebyshev2.cheb2Ord(hpSpecs);
        TransferFunction cb2 = Chebyshev2.newHighPass(hpr.getOrder(), hpSpecs.getStopBandAttenuation(), hpr.getCutoffFrequency());

        hpr = Elliptic.ellipOrd(hpSpecs);
        TransferFunction el = Elliptic.newHighPass(hpr.getOrder(), hpSpecs.getPassBandRipple(),
                hpSpecs.getStopBandAttenuation(), hpr.getCutoffFrequency());

        System.out.println();
        System.out.println("//////////////////////////////////");
        System.out.println("//");
        System.out.println("// High pass filter approximations");
        System.out.println("//");
        System.out.println("//////////////////////////////////");

        printTransferFunctions(bu, cb1, cb2, el);
    }

    public static void buildBandPassFilters(BandpassSpecs bpSpecs) {
        BandpassResults bpr = ButterWorth.buttOrd(bpSpecs);
        TransferFunction bu = ButterWorth.newBandpass(bpr.getOrder(), bpr.getLowerCutoffFrequency(),
                bpr.getUpperCutoffFrequency());

        bpr = Chebyshev1.cheb1Ord(bpSpecs);
        TransferFunction cb1 = Chebyshev1.newBandpass(bpr.getOrder(), bpSpecs.getPassBandRipple(),
                bpr.getLowerCutoffFrequency(), bpr.getUpperCutoffFrequency());

        bpr = Chebyshev2.cheb2Ord(bpSpecs);
        TransferFunction cb2 = Chebyshev2.newBandpass(bpr.getOrder(), bpSpecs.getStopBandAttenuation(),
                bpr.getLowerCutoffFrequency(), bpr.getUpperCutoffFrequency());

        bpr = Elliptic.ellipOrd(bpSpecs);
        TransferFunction el = Elliptic.newBandpass(bpr.getOrder(), bpSpecs.getPassBandRipple(),
                bpSpecs.getStopBandAttenuation(), bpr.getLowerCutoffFrequency(), bpr.getUpperCutoffFrequency());

        System.out.println();
        System.out.println("//////////////////////////////////");
        System.out.println("//");
        System.out.println("// Band pass filter approximations");
        System.out.println("//");
        System.out.println("//////////////////////////////////");

        printTransferFunctions(bu, cb1, cb2, el);
    }

    public static void buildBandStopFilters(BandStopSpecs bsSpecs) {
        BandStopResults bsr = ButterWorth.buttOrd(bsSpecs);
        TransferFunction bu = ButterWorth.newBandStop(bsr.getOrder(), bsr.getLowerCutoffFrequency(),
                bsr.getUpperCutoffFrequency());

        bsr = Chebyshev1.cheb1Ord(bsSpecs);
        TransferFunction cb1 = Chebyshev1.newBandStop(bsr.getOrder(), bsSpecs.getPassBandRipple(),
                bsr.getLowerCutoffFrequency(), bsr.getUpperCutoffFrequency());

        bsr = Chebyshev2.cheb2Ord(bsSpecs);
        TransferFunction cb2 = Chebyshev2.newBandStop(bsr.getOrder(), bsSpecs.getStopBandAttenuation(),
                bsr.getLowerCutoffFrequency(), bsr.getUpperCutoffFrequency());

        bsr = Elliptic.ellipOrd(bsSpecs);
        TransferFunction el = Elliptic.newBandStop(bsr.getOrder(), bsSpecs.getPassBandRipple(),
                bsSpecs.getStopBandAttenuation(), bsr.getLowerCutoffFrequency(), bsr.getUpperCutoffFrequency());

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
