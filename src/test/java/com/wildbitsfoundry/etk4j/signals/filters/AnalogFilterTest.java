package com.wildbitsfoundry.etk4j.signals.filters;

import static org.junit.Assert.assertArrayEquals;

import com.wildbitsfoundry.etk4j.control.ZeroPoleGain;
import org.junit.Test;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.LowPassSpecs;

import javax.swing.*;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class AnalogFilterTest {

    @Test
    public void testButterworth() {

        final LowPassSpecs lpSpecs = new LowPassSpecs();
        final HighPassSpecs hpSpecs = new HighPassSpecs();
        final BandPassSpecs bpSpecs = new BandPassSpecs();
        final BandStopSpecs bsSpecs = new BandStopSpecs();
        // TODO make sure that the filter throws if any of the input parameters are negative
        // This needs to be documented as well as a javadoc. Also that rp and rs are in db
        final double passBandRipple = -20.0 * Math.log10(1.0 / Math.sqrt(2.0));
        lpSpecs.setPassBandRipple(passBandRipple);
        lpSpecs.setStopBandAttenuation(70.0);
        lpSpecs.setPassBandFrequency(1.0 / (2.0 * Math.PI));
        lpSpecs.setStopBandFrequency(10.0 / (2.0 * Math.PI));

        hpSpecs.setPassBandRipple(passBandRipple);
        hpSpecs.setStopBandAttenuation(70.0);
        hpSpecs.setPassBandFrequency(1.0 / (2.0 * Math.PI));
        hpSpecs.setStopBandFrequency(0.1 / (2.0 * Math.PI));

        bpSpecs.setLowerPassBandFrequency(190.0);
        bpSpecs.setUpperPassBandFrequency(210.0);
        bpSpecs.setLowerStopBandFrequency(180.0);
        bpSpecs.setUpperStopBandFrequency(220.0);
        bpSpecs.setPassBandRipple(passBandRipple);
        bpSpecs.setStopBandAttenuation(20.0);

        bsSpecs.setLowerPassBandFrequency(3.6e3);
        bsSpecs.setUpperPassBandFrequency(9.1e3);
        bsSpecs.setLowerStopBandFrequency(5.45e3);
        bsSpecs.setUpperStopBandFrequency(5.90e3);
        bsSpecs.setPassBandRipple(passBandRipple);
        bsSpecs.setStopBandAttenuation(38.0);

        FilterOrderResults.OrderAndCutoffFrequency nWn = ButterWorth.buttord(lpSpecs);
        TransferFunction ba = ButterWorth.newLowPass(nWn.getOrder(), nWn.getCutoffFrequency());
        assertArrayEquals(new double[]{6.416238909177711E-4}, ba.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 0.4158919086417555, 0.08648303983684116, 0.010534665112713422,
                6.416238909177713E-4}, ba.getDenominator().getCoefficients(), 1e-12);

        nWn = ButterWorth.buttord(hpSpecs);
        ba = ButterWorth.newHighPass(nWn.getOrder(), nWn.getCutoffFrequency());
        assertArrayEquals(new double[]{1.0, 0.0, 0.0, 0.0, 0.0}, ba.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 0.4158919086417555, 0.08648303983684116, 0.010534665112713422,
                6.416238909177713E-4}, ba.getDenominator().getCoefficients(), 1e-12);

        FilterOrderResults.OrderAndCutoffFrequencies nW0W1 = ButterWorth.buttord(bpSpecs);
        ba = ButterWorth.newBandPass(nW0W1.getOrder(), nW0W1.getLowerCutoffFrequency(), nW0W1.getUpperCutoffFrequency());
        assertArrayEquals(new double[]{160000.0, 0.0, 0.0, 0.0, 0.0}, ba.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 52.26251859505506, 160965.68542494922, 6276728.483266111,
                9.66120169691095E9, 2.5044146648231787E11, 2.562589808533734E14, 3.319777843917693E15,
                2.5344958400999997E18}, ba.getDenominator().getCoefficients(), 1e-12);

        nW0W1 = ButterWorth.buttord(bsSpecs);
        ba = ButterWorth.newBandStop(nW0W1.getOrder(), nW0W1.getLowerCutoffFrequency(), nW0W1.getUpperCutoffFrequency());
        assertArrayEquals(new double[]{1.0, 0.0, 6.431000006042599E7, 0.0, 1.0339440269429976E15},
                ba.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 7540.507726234662, 9.273962844512829E7, 2.424650261648968E11,
                1.0339440269429975E15}, ba.getDenominator().getCoefficients(), 1e-12);
    }

    @Test
    public void testCheby() {

        final LowPassSpecs lpSpecs = new LowPassSpecs();
        final HighPassSpecs hpSpecs = new HighPassSpecs();
        final BandPassSpecs bpSpecs = new BandPassSpecs();
        final BandStopSpecs bsSpecs = new BandStopSpecs();

        final double passBandRipple = -20.0 * Math.log10(1.0 / Math.sqrt(2.0));
        lpSpecs.setPassBandRipple(passBandRipple);
        lpSpecs.setStopBandAttenuation(70.0);
        lpSpecs.setPassBandFrequency(1.0 / (2.0 * Math.PI));
        lpSpecs.setStopBandFrequency(10.0 / (2.0 * Math.PI));

        hpSpecs.setPassBandRipple(passBandRipple);
        hpSpecs.setStopBandAttenuation(70.0);
        hpSpecs.setPassBandFrequency(1.0 / (2.0 * Math.PI));
        hpSpecs.setStopBandFrequency(0.1 / (2.0 * Math.PI));

        bpSpecs.setLowerPassBandFrequency(190.0);
        bpSpecs.setUpperPassBandFrequency(210.0);
        bpSpecs.setLowerStopBandFrequency(180.0);
        bpSpecs.setUpperStopBandFrequency(220.0);
        bpSpecs.setPassBandRipple(passBandRipple);
        bpSpecs.setStopBandAttenuation(20.0);

        bsSpecs.setLowerPassBandFrequency(3.6e3);
        bsSpecs.setUpperPassBandFrequency(9.1e3);
        bsSpecs.setLowerStopBandFrequency(5.45e3);
        bsSpecs.setUpperStopBandFrequency(5.90e3);
        bsSpecs.setPassBandRipple(passBandRipple);
        bsSpecs.setStopBandAttenuation(38.0);

        FilterOrderResults.OrderAndCutoffFrequency nWn = Chebyshev1.cheb1ord(lpSpecs);
        TransferFunction ba = Chebyshev1.newLowPass(nWn.getOrder(), lpSpecs.getPassBandRipple(),
                nWn.getCutoffFrequency());
        assertArrayEquals(new double[]{0.0010078604510374838}, ba.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 0.09486774762192837, 0.023497666702367283, 0.0010078604510374838},
                ba.getDenominator().getCoefficients(), 1e-12);

        nWn = Chebyshev1.cheb1ord(hpSpecs);
        ba = Chebyshev1.newHighPass(nWn.getOrder(), hpSpecs.getPassBandRipple(),
                nWn.getCutoffFrequency());
        assertArrayEquals(new double[]{0.9999999999999999, 0.0, 0.0, 0.0}, ba.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 0.5905607767092684, 0.0603946838960972, 0.016125767216599744},
                ba.getDenominator().getCoefficients(), 1e-12);

        FilterOrderResults.OrderAndCutoffFrequencies nW0W1 = Chebyshev1.cheb1ord(bpSpecs);
        ba = Chebyshev1.newBandPass(nW0W1.getOrder(), bpSpecs.getPassBandRipple(), nW0W1.getLowerCutoffFrequency(),
                nW0W1.getUpperCutoffFrequency());
        assertArrayEquals(new double[]{1999.9999999999993, 0.0, 0.0, 0.0}, ba.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 11.921432759666429, 120071.06027952161, 953330.334221381,
                4.790835305152913E9, 1.8979040167716553E10, 6.3521199000000016E13}, ba.getDenominator().getCoefficients(), 1e-12);

        nW0W1 = Chebyshev1.cheb1ord(bsSpecs);
        ba = Chebyshev1.newBandStop(nW0W1.getOrder(), bsSpecs.getPassBandRipple(), nW0W1.getLowerCutoffFrequency(),
                nW0W1.getUpperCutoffFrequency());
        assertArrayEquals(new double[]{0.7071067811865476, 0.0, 4.547403706833129E7, 0.0, 7.311088304873791E14},
                ba.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 4853.027422258612, 1.0451556575522438E8, 1.5604909666054745E11,
                1.0339440236459848E15}, ba.getDenominator().getCoefficients(), 1e-12);
    }

    @Test
    public void testInverseCheby() {

        final LowPassSpecs lpSpecs = new LowPassSpecs();
        final HighPassSpecs hpSpecs = new HighPassSpecs();
        final BandPassSpecs bpSpecs = new BandPassSpecs();
        final BandStopSpecs bsSpecs = new BandStopSpecs();

        final double passBandRipple = -20.0 * Math.log10(1.0 / Math.sqrt(2.0));
        lpSpecs.setPassBandRipple(passBandRipple);
        lpSpecs.setStopBandAttenuation(70.0);
        lpSpecs.setPassBandFrequency(1.0 / (2.0 * Math.PI));
        lpSpecs.setStopBandFrequency(10.0 / (2.0 * Math.PI));

        hpSpecs.setPassBandRipple(passBandRipple);
        hpSpecs.setStopBandAttenuation(70.0);
        hpSpecs.setPassBandFrequency(1.0 / (2.0 * Math.PI));
        hpSpecs.setStopBandFrequency(0.1 / (2.0 * Math.PI));

        bpSpecs.setLowerPassBandFrequency(190.0);
        bpSpecs.setUpperPassBandFrequency(210.0);
        bpSpecs.setLowerStopBandFrequency(180.0);
        bpSpecs.setUpperStopBandFrequency(220.0);
        bpSpecs.setPassBandRipple(passBandRipple);
        bpSpecs.setStopBandAttenuation(20.0);

        bsSpecs.setLowerPassBandFrequency(3.6e3);
        bsSpecs.setUpperPassBandFrequency(9.1e3);
        bsSpecs.setLowerStopBandFrequency(5.45e3);
        bsSpecs.setUpperStopBandFrequency(5.90e3);
        bsSpecs.setPassBandRipple(passBandRipple);
        bsSpecs.setStopBandAttenuation(38.0);

        FilterOrderResults.OrderAndCutoffFrequency nWn = Chebyshev2.cheb2ord(lpSpecs);
        TransferFunction ba = Chebyshev2.newLowPass(nWn.getOrder(), lpSpecs.getStopBandAttenuation(),
                nWn.getCutoffFrequency());
        assertArrayEquals(new double[]{0.0014001979322511749, 0.0, 0.004066909232107251}, ba.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 0.3187766294691388, 0.050808289470727575, 0.004066909232107251},
                ba.getDenominator().getCoefficients(), 1e-12);

        nWn = Chebyshev2.cheb2ord(hpSpecs);
        ba = Chebyshev2.newHighPass(nWn.getOrder(), hpSpecs.getStopBandAttenuation(),
                nWn.getCutoffFrequency());
        assertArrayEquals(new double[]{0.9999999999999997, 0.0, 2.2090496592679388E-4, 0.0}, ba.getNumerator().getCoefficients(),
                1e-12);
        assertArrayEquals(new double[]{1.0, 0.3164538310429205, 0.05029241855679742, 0.003996283686869137},
                ba.getDenominator().getCoefficients(), 1e-12);

        FilterOrderResults.OrderAndCutoffFrequencies nW0W1 = Chebyshev2.cheb2ord(bpSpecs);
        ba = Chebyshev2.newBandPass(nW0W1.getOrder(), bpSpecs.getStopBandAttenuation(), nW0W1.getLowerCutoffFrequency(),
                nW0W1.getUpperCutoffFrequency());
        assertArrayEquals(new double[]{9.277258525655208, 0.0, 752036.1337575477, 0.0, 1.4769488345428343E10, 0.0},
                ba.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 43.24250812255668, 120591.92349148876, 3462463.051590285,
                4.811617747310401E9, 6.884250535619144E10, 6.352119899999999E13}, ba.getDenominator().getCoefficients(), 1e-12);

        nW0W1 = Chebyshev2.cheb2ord(bsSpecs);
        ba = Chebyshev2.newBandStop(nW0W1.getOrder(), bsSpecs.getStopBandAttenuation(), nW0W1.getLowerCutoffFrequency(),
                nW0W1.getUpperCutoffFrequency());
        assertArrayEquals(new double[]{1.0, 0.0, 6.466348566357853E7, 0.0, 1.0339440236459846E15},
                ba.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 7446.459823865073, 9.238836761779684E7, 2.3944091547959982E11,
                1.0339440236459841E15}, ba.getDenominator().getCoefficients(), 1e-12);
    }

    @Test
    public void testElliptic() {

        final LowPassSpecs lpSpecs = new LowPassSpecs();
        final HighPassSpecs hpSpecs = new HighPassSpecs();
        final BandPassSpecs bpSpecs = new BandPassSpecs();
        final BandStopSpecs bsSpecs = new BandStopSpecs();

        final double passBandRipple = -20.0 * Math.log10(1.0 / Math.sqrt(2.0));
        lpSpecs.setPassBandRipple(passBandRipple);
        lpSpecs.setStopBandAttenuation(70.0);
        lpSpecs.setPassBandFrequency(1.0 / (2.0 * Math.PI));
        lpSpecs.setStopBandFrequency(10.0 / (2.0 * Math.PI));

        hpSpecs.setPassBandRipple(passBandRipple);
        hpSpecs.setStopBandAttenuation(70.0);
        hpSpecs.setPassBandFrequency(1.0 / (2.0 * Math.PI));
        hpSpecs.setStopBandFrequency(0.1 / (2.0 * Math.PI));

        bpSpecs.setLowerPassBandFrequency(190.0);
        bpSpecs.setUpperPassBandFrequency(210.0);
        bpSpecs.setLowerStopBandFrequency(180.0);
        bpSpecs.setUpperStopBandFrequency(220.0);
        bpSpecs.setPassBandRipple(passBandRipple);
        bpSpecs.setStopBandAttenuation(20.0);

        bsSpecs.setLowerPassBandFrequency(3.6e3);
        bsSpecs.setUpperPassBandFrequency(9.1e3);
        bsSpecs.setLowerStopBandFrequency(5.45e3);
        bsSpecs.setUpperStopBandFrequency(5.90e3);
        bsSpecs.setPassBandRipple(passBandRipple);
        bsSpecs.setStopBandAttenuation(38.0);

        FilterOrderResults.OrderAndCutoffFrequency nWn = Elliptic.ellipord(lpSpecs);
        TransferFunction ba = Elliptic.newLowPass(nWn.getOrder(), lpSpecs.getPassBandRipple(),
                lpSpecs.getStopBandAttenuation(), nWn.getCutoffFrequency());
        assertArrayEquals(new double[]{8.794918109662971E-4, 0.0, 0.0010190301217685481}, ba.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 0.09475123472716065, 0.02355602777927699, 0.0010190301217685483},
                ba.getDenominator().getCoefficients(), 1e-12);

        nWn = Elliptic.ellipord(hpSpecs);
        ba = Elliptic.newHighPass(nWn.getOrder(), lpSpecs.getPassBandRipple(),
                lpSpecs.getStopBandAttenuation(), nWn.getCutoffFrequency());
        assertArrayEquals(new double[]{1.0, 0.0, 5.5376474721194E-4, 0.0}, ba.getNumerator().getCoefficients(),
                1e-12);
        assertArrayEquals(new double[]{1.0, 0.58553828918372, 0.059659331550860784, 0.015949011391381742},
                ba.getDenominator().getCoefficients(), 1e-12);

        FilterOrderResults.OrderAndCutoffFrequencies nW0W1 = Elliptic.ellipord(bpSpecs);
        ba = Elliptic.newBandPass(nW0W1.getOrder(), bpSpecs.getPassBandRipple(), bpSpecs.getStopBandAttenuation(),
                nW0W1.getLowerCutoffFrequency(), nW0W1.getUpperCutoffFrequency());
        assertArrayEquals(new double[]{0.10000000000000003, 0.0, 8198.997487421326, 0.0, 1.5920100000000003E8},
                ba.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 11.824526997102488, 80109.70921683687, 471798.6271843893,
                1.5920099999999998E9}, ba.getDenominator().getCoefficients(), 1e-12);

        nW0W1 = Elliptic.ellipord(bsSpecs);
        ba = Elliptic.newBandStop(nW0W1.getOrder(), bsSpecs.getPassBandRipple(), bsSpecs.getStopBandAttenuation(),
                nW0W1.getLowerCutoffFrequency(), nW0W1.getUpperCutoffFrequency());
        assertArrayEquals(new double[]{0.7071067811865475, 0.0, 4.597394172926564E7, 0.0, 7.311088423294805E14},
                ba.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 4749.564283857035, 1.0401880991009936E8, 1.527222406842788E11,
                1.0339440403932445E15}, ba.getDenominator().getCoefficients(), 1e-12);
    }

    @Test
    public void testBessel() {

        TransferFunction ba = Bessel.newLowPass(2, 1000);
        assertArrayEquals(new double[]{1000000.0}, ba.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 1732.0508075688772, 999999.9999999995}, ba.getDenominator().getCoefficients(), 1e-12);

        ba = Bessel.newHighPass(2, 1000);
        assertArrayEquals(new double[]{1.0000000000000004, 0.0, 0.0}, ba.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 1732.0508075688779, 1000000.0000000005}, ba.getDenominator().getCoefficients(), 1e-12);

        ba = Bessel.newBandPass(4, 190e6, 210e6);
        assertArrayEquals(new double[]{1.6E29, 0.0, 0.0, 0.0, 0.0}, ba.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 6.2478798738405116E7, 1.61356620131307328E17, 7.504320895970641E24,
                9.692398286478327E33, 2.9942240374922867E41, 2.5688135281524276E50, 3.968728207943182E57,
                2.534495840100001E66}, ba.getDenominator().getCoefficients(), 1e-12);

        ba = Bessel.newBandStop(4, 180e6, 220e6);
        assertArrayEquals(new double[]{1.0000000000000004, 0.0, 1.584E17, 0.0, 9.408960000000003E33, 0.0,
                2.4839654400000013E50, 0.0, 2.459125785600001E66},
                ba.getNumerator().getCoefficients(), 1e-12);
        assertArrayEquals(new double[]{1.0, 1.280434349177472E8, 1.6542648052522944E17, 1.5411492224191265E25,
                9.968017257598169E33, 6.1029509207797396E41, 2.5941518970044378E50, 7.95138667886433E57,
                2.4591257855999993E66}, ba.getDenominator().getCoefficients(), 1e-12);
    }
}
