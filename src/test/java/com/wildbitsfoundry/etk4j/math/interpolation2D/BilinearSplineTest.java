package com.wildbitsfoundry.etk4j.math.interpolation2D;

import static com.wildbitsfoundry.etk4j.math.interpolation2D.Spline2D.newBilinearSpline;
import static org.junit.Assert.assertEquals;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

public class BilinearSplineTest {
	static double[] x;
	static double[] y;
	static double[][] z;
	static double[][] xiyi;

	@Rule
	public final ExpectedException exception = ExpectedException.none();

	@BeforeClass
	public static void setArrays() {
		x = new double[] { 1, 2, 3, 4, 5, 6, 7, 8 };
		y = new double[] { 1, 2, 3, 4, 5, 6, 7, 8 };

		z = new double[][] { { 1, 4, 9, 16, 25, 36, 49, 64 }, { 4, 16, 36, 64, 100, 144, 196, 256 },
				{ 9, 36, 81, 144, 225, 324, 441, 576 }, { 16, 64, 144, 256, 400, 576, 784, 1024 },
				{ 25, 100, 225, 400, 625, 900, 1225, 1600 }, { 36, 144, 324, 576, 900, 1296, 1764, 2304 },
				{ 49, 196, 441, 784, 1225, 1764, 2401, 3136 }, { 64, 256, 576, 1024, 1600, 2304, 3136, 4096 } };
		xiyi = new double[][] { { 1.5, 1.5 }, { 2.5, 2.5 }, { 3.5, 3.5 }, { 4.5, 4.5 }, { 5.5, 5.5 }, { 6.5, 6.5 },
				{ 7.5, 7.5 } };
	}

	@Test
	public void testBilinearSplineInterpolation() {
		Spline2D sp = newBilinearSpline(x, y, z);

		double yi = sp.evaluateAt(xiyi[0][0], xiyi[0][1]);
		assertEquals(6.2500, yi, 1e-12);
		
		yi = sp.evaluateAt(xiyi[1][0], xiyi[1][1]);
		assertEquals(42.2500, yi, 1e-12);
		
		yi = sp.evaluateAt(xiyi[2][0], xiyi[2][1]);
		assertEquals(156.2500, yi, 1e-12);
		
		yi = sp.evaluateAt(xiyi[3][0], xiyi[3][1]);
		assertEquals(420.2500, yi, 1e-12);
		
		yi = sp.evaluateAt(xiyi[4][0], xiyi[4][1]);
		assertEquals(930.2500, yi, 1e-12);
		
		yi = sp.evaluateAt(xiyi[5][0], xiyi[5][1]);
		assertEquals(1806.2500, yi, 1e-12);
		
		yi = sp.evaluateAt(xiyi[6][0], xiyi[6][1]);
		assertEquals(3192.2500, yi, 1e-12);
	}
}