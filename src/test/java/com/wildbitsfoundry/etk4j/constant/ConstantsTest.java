package com.wildbitsfoundry.etk4j.constant;

import static org.junit.Assert.assertEquals;

import org.junit.Test;
import static com.wildbitsfoundry.etk4j.constant.ConstantsETK.*;

public class ConstantsTest {

	@Test
	public void testAll() {
		assertEquals(2.220446049250313E-16, DOUBLE_EPS, 1e-16);
		assertEquals((float) 1.1920929E-7, FLOAT_EPS, 1e-8);
		assertEquals(6.283185307179586, TWO_PI, 1e-16);
	}
}
