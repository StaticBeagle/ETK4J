package com.wildbitsfoundry.etk4j.math.calculus;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import org.junit.Test;

public final class IntegralsTest {
	
	double[] a = {1, 2, 3, 4 ,5 ,6};
	
	@Test
	public void testCummulativeTrapz() {
		double[] trapz = Integrals.cummulativeTrapz(a);
		assertArrayEquals(new double[] {0.0, 1.5, 4.0, 7.5, 12.0, 17.5}, trapz, 1e-12);
	}
	
	@Test
	public void testTrapz() {
		double trapz = Integrals.trapz(a);
		assertEquals(17.5, trapz, 1e-12);
	}	
}
