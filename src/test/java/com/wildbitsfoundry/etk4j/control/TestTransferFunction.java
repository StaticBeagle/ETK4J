package com.wildbitsfoundry.etk4j.control;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.junit.Test;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.util.NumArrays;

public class TestTransferFunction {

	@Test
	public void testConstructors() {
		TransferFunction tf = new TransferFunction(new double[] { 10.0 }, new double[] { 1.0, 1.0 });
		assertEquals(10.0, tf.getMagnitudeAt(0.0), 1e-12);

	}

	@Test
	public void testMargins() {
		double[] numerator = { 1.0 };
		double[] denominator = { 1.0, 2.0, 1.0, 0.0 };

		TransferFunction tf = new TransferFunction(numerator, denominator);
		Margins margins = tf.getMargins();

		assertEquals(0.6823278038280193, margins.getGainCrossOverFrequency(), 1e-12);
		assertEquals(21.386389751875015, margins.getPhaseMargin(), 1e-12);
		assertEquals(0.9999999999999998, margins.getPhaseCrossOverFrequency(), 1e-12);
		assertEquals(6.02059991327962, margins.getGainMargin(), 1e-12);
		assertEquals(
				"Margins [GainMargin=6.02059991327962, PhaseMargin=21.386389751875015,"
						+ " GainCrossOverFrequency=0.6823278038280193, PhaseCrossOverFrequency=0.9999999999999998]",
				margins.toString());
	}

	@Test
	public void testPolesAndZeros() {
		double[] numerator = { 1.0 };
		double[] denominator = { 1.0, 2.0, 1.0, 0.0 };

		TransferFunction tf = new TransferFunction(numerator, denominator);

		Complex[] zeros = tf.getZeros();
		Complex[] poles = tf.getPoles();

		assertArrayEquals(new Complex[] {}, zeros);
		assertArrayEquals(new Complex[] { Complex.newComplex(-1.0000000209081399, 0.0), Complex.newComplex(-0.9999999790918601, 0.0),
				Complex.newComplex(-4.1910912494273124E-17, 0.0) }, poles);

		zeros = new Complex[] { Complex.newComplex(-1.0, 0.0) };
		poles = new Complex[] { Complex.newComplex(-5.0, 0.0), Complex.newComplex(-1.0, 0.0), Complex.newComplex(-2.0, 0.0) };

		tf = new TransferFunction(zeros, poles);

		numerator = tf.getNumerator().getCoefficients();
		denominator = tf.getDenominator().getCoefficients();

		assertArrayEquals(new double[] { 10.0, 10.0 }, numerator, 1e-12);
		assertArrayEquals(new double[] { 1.0, 8.0, 17.0, 10.0 }, denominator, 1e-12);
	}

	@Test
	public void testGettersAndEvaluation() {
		Complex[] poles = new Complex[] { Complex.newComplex(-1.0, 0.0), Complex.newComplex(-1.0, 0.0), Complex.newComplex(-1.0, 0.0) };

		TransferFunction tf = new TransferFunction(10.0, poles);

		double phase = tf.getPhaseAt(100.0);
		assertEquals(-268.2811839069496, phase, 1e-12);

		double[] frequencies = NumArrays.logspace(-3, 3, 10);

		double[] magnitudeResponse = { 9.999985000018752, 9.999676843499255, 9.993041654128266, 9.851853368415734,
				7.462732134984385, 0.7462732134984399, 0.009851853368415734, 9.993041654128302E-5, 9.999676843499274E-7,
				9.999985000018753E-9 };

		double[] phaseResponse = { -0.17188728124350178, -0.7978246216992994, -3.7026276509801344, -17.13177941249893,
				-74.69637093434646, -195.30362906565347, -252.86822058750107, -266.29737234901984, -269.2021753783007,
				-269.8281127187565 };
		
		assertArrayEquals(magnitudeResponse, tf.getMagnitudeAt(frequencies), 1e-12);
		assertArrayEquals(phaseResponse, tf.getPhaseAt(frequencies), 1e-12);
	}

}
