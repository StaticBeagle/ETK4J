package com.wildbitsfoundry.etk4j.signals.filters;

import static org.junit.Assert.assertArrayEquals;

import org.junit.BeforeClass;
import org.junit.Test;

import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.BandStopSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.HighPassSpecs;
import com.wildbitsfoundry.etk4j.signals.filters.FilterSpecs.LowPassSpecs;

public class AnalogFilterTest {
    static LowPassSpecs lpSpecs = new LowPassSpecs();
    static HighPassSpecs hpSpecs = new HighPassSpecs();
    static BandPassSpecs bpSpecs = new BandPassSpecs();
    static BandStopSpecs bsSpecs = new BandStopSpecs();

    @BeforeClass
    public static void setUpClass() {
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
        bpSpecs.setLowerStopBandAttenuation(20);
        bpSpecs.setUpperStopBandAttenuation(20);

        bsSpecs.setLowerPassBandFrequency(3.6e3);
        bsSpecs.setUpperPassBandFrequency(9.1e3);
        bsSpecs.setLowerStopBandFrequency(5.45e3);
        bsSpecs.setUpperStopBandFrequency(5.90e3);
        bsSpecs.setPassBandRipple(passBandRipple);
        bsSpecs.setStopBandAttenuation(38.0);
    }

    @Test
    public void testButterworth() {

        AnalogFilter lp = AnalogFilter.newLowPass(lpSpecs, ApproximationType.BUTTERWORTH);
        assertArrayEquals(new double[]{6.416238909177711E-4}, lp.getNumerator(), 1e-12);
        assertArrayEquals(new double[]{1.0, 0.4158919086417556, 0.08648303983684116, 0.01053466511271342,
						6.416238909177711E-4}, lp.getDenominator(), 1e-12);

        AnalogFilter hp = AnalogFilter.newHighPass(hpSpecs, ApproximationType.BUTTERWORTH);
        assertArrayEquals(new double[]{1.0, 0.0, 0.0, 0.0, 0.0}, hp.getNumerator(), 1e-12);
        assertArrayEquals(new double[]{1.0, 0.41589190864175557, 0.08648303983684116, 0.01053466511271342,
						6.416238909177711E-4}, hp.getDenominator(), 1e-12);

        AnalogFilter bp = AnalogFilter.newBandPass(bpSpecs, ApproximationType.BUTTERWORTH);
        assertArrayEquals(new double[]{160000.0, 0.0, 0.0, 0.0, 0.0}, bp.getNumerator(), 1e-12);
        assertArrayEquals(new double[]{1.0, 52.26251859505506, 160965.68542494925, 6276728.483266112,
				9.66120169691095E9, 2.504414664823179E11, 2.5625898085337344E14, 3.319777843917693E15,
				2.5344958400999997E18}, bp.getDenominator(), 1e-12);

        AnalogFilter bs = AnalogFilter.newBandStop(bsSpecs, ApproximationType.BUTTERWORTH);
        assertArrayEquals(new double[]{1.0, 0.0, 6.431000004679801E7, 0.0, 1.03394402650479E15},
				bs.getNumerator(), 1e-12);
        assertArrayEquals(new double[]{1.0, 7540.508141142473, 9.273963156011595E7,
				2.4246503945487656E11, 1.0339440265047899E15}, bs.getDenominator(), 1e-12);
    }

    @Test
    public void testCheby() {

        AnalogFilter lp = AnalogFilter.newLowPass(lpSpecs, ApproximationType.CHEBYSHEV);
        assertArrayEquals(new double[]{0.0010078604510374838}, lp.getNumerator(), 1e-12);
        assertArrayEquals(new double[]{1.0, 0.09486774762192837, 0.02349766670236728, 0.0010078604510374838},
                lp.getDenominator(), 1e-12);

        AnalogFilter hp = AnalogFilter.newHighPass(hpSpecs, ApproximationType.CHEBYSHEV);
        assertArrayEquals(new double[]{1.0, 0.0, 0.0, 0.0}, hp.getNumerator(), 1e-12);
        assertArrayEquals(new double[]{1.0, 0.5905607767092685, 0.06039468389609723, 0.016125767216599755},
                hp.getDenominator(), 1e-12);

        AnalogFilter bp = AnalogFilter.newBandPass(bpSpecs, ApproximationType.CHEBYSHEV);
        assertArrayEquals(new double[]{1999.9999999999993, 0.0, 0.0, 0.0}, bp.getNumerator(), 1e-12);
        assertArrayEquals(new double[]{1.0, 11.92143275966643, 120071.06027952163, 953330.3342213811,
                4.790835305152913E9, 1.8979040167716553E10, 6.3521199E13}, bp.getDenominator(), 1e-12);

        AnalogFilter bs = AnalogFilter.newBandStop(bsSpecs, ApproximationType.CHEBYSHEV);
        assertArrayEquals(new double[]{0.7071067811865475, 0.0, 4.5474037131198056E7, 0.0, 7.311088325088602E14},
                bs.getNumerator(), 1e-12);
        assertArrayEquals(new double[]{1.0, 4853.0277036270545, 1.0451557050620179E8,
                1.5604910592368393E11, 1.0339440265047898E15}, bs.getDenominator(), 1e-12);
    }

    @Test
    public void testInverseCheby() {

        AnalogFilter lp = AnalogFilter.newLowPass(lpSpecs, ApproximationType.INVERSE_CHEBYSHEV);
        assertArrayEquals(new double[]{0.0014001979322511755, 0.0, 0.004066909232107252}, lp.getNumerator(), 1e-12);
        assertArrayEquals(new double[]{1.0, 0.3187766294691389, 0.05080828947072759, 0.0040669092321072514},
                lp.getDenominator(), 1e-12);

        AnalogFilter hp = AnalogFilter.newHighPass(hpSpecs, ApproximationType.INVERSE_CHEBYSHEV);
        assertArrayEquals(new double[]{1.0000000000000002, 0.0, 2.2090496592679413E-4, 0.0}, hp.getNumerator(),
                1e-12);
        assertArrayEquals(new double[]{1.0, 0.3164538310429205, 0.05029241855679743, 0.003996283686869137},
                hp.getDenominator(), 1e-12);

        AnalogFilter bp = AnalogFilter.newBandPass(bpSpecs, ApproximationType.INVERSE_CHEBYSHEV);
        assertArrayEquals(new double[]{9.277258525655194, 0.0, 752036.1337575465, 0.0, 1.4769488345428324E10, 0.0},
                bp.getNumerator(), 1e-12);
        assertArrayEquals(new double[]{1.0, 43.2425081225566, 120591.92349148876, 3462463.0515902787,
                4.811617747310402E9, 6.884250535619133E10, 6.3521199E13}, bp.getDenominator(), 1e-12);

        AnalogFilter bs = AnalogFilter.newBandStop(bsSpecs, ApproximationType.INVERSE_CHEBYSHEV);
        assertArrayEquals(new double[]{1.0, 0.0, 6.4663485793474294E7, 0.0, 1.0339440265047899E15},
                bs.getNumerator(), 1e-12);
        assertArrayEquals(new double[]{1.0, 7446.460255595333, 9.238837096255475E7,
                2.3944092969290768E11, 1.0339440265047899E15}, bs.getDenominator(), 1e-12);
    }

    @Test
    public void testElliptic() {

        AnalogFilter lp = AnalogFilter.newLowPass(lpSpecs, ApproximationType.ELLIPTIC);
        assertArrayEquals(new double[]{8.794918109662971E-4, 0.0, 0.0010190301217685481}, lp.getNumerator(), 1e-12);
        assertArrayEquals(new double[]{1.0, 0.09475123472716064, 0.023556027779276987, 0.0010190301217685481},
                lp.getDenominator(), 1e-12);

        AnalogFilter hp = AnalogFilter.newHighPass(hpSpecs, ApproximationType.ELLIPTIC);
        assertArrayEquals(new double[]{1.0, 0.0, 5.537647472119399E-4, 0.0}, hp.getNumerator(),
                1e-12);
        assertArrayEquals(new double[]{1.0, 0.5855382891837199, 0.059659331550860784, 0.015949011391381742},
                hp.getDenominator(), 1e-12);

        AnalogFilter bp = AnalogFilter.newBandPass(bpSpecs, ApproximationType.ELLIPTIC);
        assertArrayEquals(new double[]{0.10000000000000003, 0.0, 8198.997487421328, 0.0, 1.5920100000000006E8},
                bp.getNumerator(), 1e-12);
        assertArrayEquals(new double[]{1.0, 11.824526997102488, 80109.70921683687, 471798.62718438933, 1.59201E9},
                bp.getDenominator(), 1e-12);

        AnalogFilter bs = AnalogFilter.newBandStop(bsSpecs, ApproximationType.ELLIPTIC);
        assertArrayEquals(new double[]{0.7071067811865474, 0.0, 4.597394146825323E7, 0.0, 7.311088325088601E14},
                bs.getNumerator(), 1e-12);
        assertArrayEquals(new double[]{1.0, 4749.564494790546, 1.0401881300520346E8,
				1.5272224644112506E11, 1.0339440265047899E15}, bs.getDenominator(), 1e-12);
    }
}
