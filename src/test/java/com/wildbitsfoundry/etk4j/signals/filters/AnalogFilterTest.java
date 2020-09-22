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
		lpSpecs.setPassBandRipple(-20.0 * Math.log10(1.0 / Math.sqrt(2.0)));
		lpSpecs.setStopBandAttenuation(70.0);
		lpSpecs.setPassBandFrequency(1.0 / (2.0 * Math.PI));
		lpSpecs.setStopBandFrequency(10.0 / (2.0 * Math.PI));

		hpSpecs.setPassBandRipple(-20.0 * Math.log10(1.0 / Math.sqrt(2.0)));
		hpSpecs.setStopBandAttenuation(70.0);
		hpSpecs.setPassBandFrequency(1.0 / (2.0 * Math.PI));
		hpSpecs.setStopBandFrequency(0.1 / (2.0 * Math.PI));

		bpSpecs.setLowerPassBandFrequency(190.0);
		bpSpecs.setUpperPassBandFrequency(210.0);
		bpSpecs.setLowerStopBandFrequency(180.0);
		bpSpecs.setUpperStopBandFrequency(220.0);
		bpSpecs.setPassBandRipple(-20.0 * Math.log10(1.0 / Math.sqrt(2.0)));
		bpSpecs.setLowerStopBandAttenuation(20);
		bpSpecs.setUpperStopBandAttenuation(20);

		bsSpecs.setLowerPassBandFrequency(3.6e3);
		bsSpecs.setUpperPassBandFrequency(9.1e3);
		bsSpecs.setLowerStopBandFrequency(5.45e3);
		bsSpecs.setUpperStopBandFrequency(5.90e3);
		bsSpecs.setPassBandRipple(-20.0 * Math.log10(1.0 / Math.sqrt(2.0)));
		bsSpecs.setStopBandAttenuation(38.0);
	}

	@Test
	public void testButterworth() {

		AnalogFilter lp = new Butterworth(lpSpecs);
		assertArrayEquals(lp.getNumerator(), new double[] { 1.0 }, 1e-12);
		assertArrayEquals(lp.getDenominator(),
				new double[] { 1.0, 2.613125929752753, 3.414213562373095, 2.613125929752753, 1.0 }, 1e-12);

		AnalogFilter hp = new Butterworth(hpSpecs);
		assertArrayEquals(hp.getNumerator(), new double[] { 1.0, 0.0, 0.0, 0.0, 0.0 }, 1e-12);
		assertArrayEquals(hp.getDenominator(),
				new double[] { 1.0, 2.613125929752753, 3.414213562373095, 2.613125929752753, 1.0 }, 1e-12);

		AnalogFilter bp = new Butterworth(bpSpecs);
		assertArrayEquals(bp.getNumerator(), new double[] { 2.4936727304704642E8, 0.0, 0.0, 0.0, 0.0 }, 1e-12);
		assertArrayEquals(bp.getDenominator(),
				new double[] { 1.0, 328.37508895264995, 6354670.549177777, 1.5569438399941423E9, 1.5057422009476117E13,
						2.452480596162642E15, 1.5767335356739457E19, 1.2834179250826303E21, 6.156452451556017E24 },
				1e-12);

		AnalogFilter bs = new Butterworth(bsSpecs);
		assertArrayEquals(bs.getNumerator(),
				new double[] { 0.9999999999999993, 0.0, 2.5866259214374967E9, 0.0, 1.67265841436309606E18 }, 1e-12);
		assertArrayEquals(bs.getDenominator(), new double[] { 0.9999999999999997, 48871.71231974204,
				3.7808480539693108E9, 6.320641895564055E13, 1.67265841436309658E18 }, 1e-12);
	}

	@Test
	public void testCheby() {

		AnalogFilter lp = new Chebyshev(lpSpecs);
		assertArrayEquals(lp.getNumerator(), new double[] { 0.06249999999999996 }, 1e-12);
		assertArrayEquals(lp.getDenominator(),
				new double[] { 1.0, 0.5960716379833215, 0.9276506988040598, 0.24999999999999992 }, 1e-12);

		AnalogFilter hp = new Chebyshev(hpSpecs);
		assertArrayEquals(hp.getNumerator(), new double[] { 0.06249999999999996, 0.0, 0.0, 0.0 }, 1e-12);
		assertArrayEquals(hp.getDenominator(),
				new double[] { 0.24999999999999992, 0.9276506988040598, 0.5960716379833215, 1.0 }, 1e-12);

		AnalogFilter bp = new Chebyshev(bpSpecs);
		assertArrayEquals(bp.getNumerator(), new double[] { 124025.10672119926, 0.0, 0.0, 0.0 }, 1e-12);
		assertArrayEquals(bp.getDenominator(), new double[] { 1.0, 74.90457115606551, 4740215.459912929,
				2.364737928847268E8, 7.466734597896846E12, 1.8585471646885812E14, 3.9083900340189711E18 }, 1e-12);

		AnalogFilter bs = new Chebyshev(bsSpecs);
		assertArrayEquals(bs.getNumerator(),
				new double[] { 0.35355339059327345, 0.0, 9.145103647206777E8, 0.0, 5.9137405370244147E17 }, 1e-12);
		assertArrayEquals(bs.getDenominator(), new double[] { 0.7071067811865472, 22241.020745028698,
				3.0232428619731693E9, 2.8764600389160188E13, 1.18274810740488346E18 }, 1e-12);
	}

	@Test
	public void testInverseCheby() {

		AnalogFilter lp = new InverseChebyshev(lpSpecs);
		assertArrayEquals(lp.getNumerator(), new double[] { 9.000000900000112E-9, 0.0, 1.2000001200000148E-6 }, 1e-12);
		assertArrayEquals(lp.getDenominator(), new double[] { 0.0010000000000000002, 0.0021598237816342143,
				0.002332374383851859, 0.0012649111273129112 }, 1e-12);

		AnalogFilter hp = new InverseChebyshev(hpSpecs);
		assertArrayEquals(hp.getNumerator(), new double[] { 1.2000001200000148E-6, 0.0, 9.000000900000109E-9, 0.0 },
				1e-12);
		assertArrayEquals(hp.getDenominator(), new double[] { 0.0012649111273129112, 0.0023323743838518586,
				0.002159823781634214, 9.999999999999998E-4 }, 1e-12);

		AnalogFilter bp = new InverseChebyshev(bpSpecs);
		assertArrayEquals(bp.getNumerator(),
				new double[] { 11.423973285781074, 0.0, 3.6230364508859456E7, 0.0, 2.8345390450910035E13, 0.0 }, 1e-12);
		assertArrayEquals(bp.getDenominator(), new double[] { 1.0, 176.6057087769558, 4740443.586016902,
				5.571724505110302E8, 7.46709393959505E12, 4.381976083026403E14, 3.9083900340189711E18 }, 1e-12);

		AnalogFilter bs = new InverseChebyshev(bsSpecs);
		assertArrayEquals(bs.getNumerator(),
				new double[] { 3.1697863849222295E-4, 0.0, 1009176.6156795565, 0.0, 5.30196986847375E14 }, 1e-12);
		assertArrayEquals(bs.getDenominator(), new double[] { 0.02517850823588335, 7705.86645274874,
				1.2593495145978765E9, 9.966096956907758E12, 4.2115043661860824E16 }, 1e-12);
	}

	@Test
	public void testElliptic() {

		AnalogFilter lp = new Elliptic(lpSpecs);
		assertArrayEquals(lp.getNumerator(), new double[] { 3.053678679030203E-5, 0.0, 0.0013968130861546817 }, 1e-12);
		assertArrayEquals(lp.getDenominator(),
				new double[] { 1.0, 0.59533956587482, 0.9299547017701413, 0.25277063920892173 }, 1e-12);

		AnalogFilter hp = new Elliptic(hpSpecs);
		assertArrayEquals(hp.getNumerator(), new double[] { 0.0013968130861546817, 0.0, 3.053678679030203E-5, 0.0 },
				1e-12);
		assertArrayEquals(hp.getDenominator(),
				new double[] { 0.25277063920892173, 0.9299547017701413, 0.59533956587482, 1.0 }, 1e-12);

		AnalogFilter bp = new Elliptic(bpSpecs);
		assertArrayEquals(bp.getNumerator(),
				new double[] { 0.014142135623730961, 0.0, 45775.752030315074, 0.0, 3.5089749077347435E10 }, 1e-12);
		assertArrayEquals(bp.getDenominator(), new double[] { 1.0, 74.29569429254273, 3162604.5546259275,
				1.1702975017491841E8, 2.4812199522726753E12 }, 1e-12);

		AnalogFilter bs = new Elliptic(bsSpecs);
		assertArrayEquals(bs.getNumerator(),
				new double[] { 0.00901331036955213, 0.0, 2.3581732495384E7, 0.0, 1.5076189430897526E16 }, 1e-12);
		assertArrayEquals(bs.getDenominator(), new double[] { 0.7159526914868412, 22039.161022784196,
				3.046123922854621E9, 2.8503532594134285E13, 1.19754429370137139E18 }, 1e-12);
	}
}
