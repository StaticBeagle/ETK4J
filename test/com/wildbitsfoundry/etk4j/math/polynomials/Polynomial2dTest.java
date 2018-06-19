package com.wildbitsfoundry.etk4j.math.polynomials;

import static com.wildbitsfoundry.etk4j.math.polynomials.Polynomial2d.polyFit2d;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import java.util.Arrays;

import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import com.wildbitsfoundry.etk4j.util.NumArrays;

public class Polynomial2dTest {
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
	public void testPolynomial2d() {
		Polynomial2d poly = polyFit2d(x, y, z, 2, 2);

		double yi = poly.evaluateAt(xiyi[0][0], xiyi[0][1]);
		assertEquals(5.0625, yi, 1e-12);

		yi = poly.evaluateAt(xiyi[1][0], xiyi[1][1]);
		assertEquals(39.0625, yi, 1e-12);

		yi = poly.evaluateAt(xiyi[2][0], xiyi[2][1]);
		assertEquals(150.0625, yi, 1e-12);

		yi = poly.evaluateAt(xiyi[3][0], xiyi[3][1]);
		assertEquals(410.0625, yi, 1e-12);

		yi = poly.evaluateAt(xiyi[4][0], xiyi[4][1]);
		assertEquals(915.0625, yi, 1e-12);

		yi = poly.evaluateAt(xiyi[5][0], xiyi[5][1]);
		assertEquals(1785.0625, yi, 1e-12);

		yi = poly.evaluateAt(xiyi[6][0], xiyi[6][1]);
		assertEquals(3164.0625, yi, 1e-12);

		double[] yyi = NumArrays.addElementWise(x, 0.5);
		double[] result = poly.evaluateAt(x, y);
	}

	@Test
	public void testPolynomial2dFlatPolyFit() {
		int length = z.length;
		double[] xx = NumArrays.repeat(x, length);
		double[] yy = NumArrays.repeat(y, length);
		Arrays.sort(yy);

		double[] zz = NumArrays.flatten(z);
		Polynomial2d poly = polyFit2d(xx, yy, zz, 2, 2);

		double yi = poly.evaluateAt(xiyi[0][0], xiyi[0][1]);
		assertEquals(5.0625, yi, 1e-12);

		yi = poly.evaluateAt(xiyi[1][0], xiyi[1][1]);
		assertEquals(39.0625, yi, 1e-12);

		yi = poly.evaluateAt(xiyi[2][0], xiyi[2][1]);
		assertEquals(150.0625, yi, 1e-12);

		yi = poly.evaluateAt(xiyi[3][0], xiyi[3][1]);
		assertEquals(410.0625, yi, 1e-12);

		yi = poly.evaluateAt(xiyi[4][0], xiyi[4][1]);
		assertEquals(915.0625, yi, 1e-12);

		yi = poly.evaluateAt(xiyi[5][0], xiyi[5][1]);
		assertEquals(1785.0625, yi, 1e-12);

		yi = poly.evaluateAt(xiyi[6][0], xiyi[6][1]);
		assertEquals(3164.0625, yi, 1e-12);
	}

	@Test
	public void testPolynomial2dArrayEval() {
		Polynomial2d poly = polyFit2d(x, y, z, 2, 2);

		double[] yi = new double[xiyi.length];
		for (int i = 0; i < xiyi.length; ++i) {
			yi[i] = xiyi[i][0];
		}
		assertArrayEquals(new double[] { 5.0625, 39.0625, 150.0625, 410.0625, 915.0625, 1785.0625, 3164.0625 },
				poly.evaluateAt(yi, yi), 1e-12);
	}

	@Test
	public void testToString() {
		Polynomial2d poly = polyFit2d(x, y, z, 2, 2);

		assertEquals("x^2 * y^2", poly.toString());

		poly = new Polynomial2d(new double[] { 1, 1, 1, 1 }, 1, 1);
		assertEquals("x * y + y + x + 1.000", poly.toString());
		
		poly = new Polynomial2d(new double[] { 2, 1, 1, 2 }, 1, 1);
		assertEquals("2.000 * x * y + y + x + 2.000", poly.toString());
		
		poly = new Polynomial2d(new double[] { 1, 1, 1, 0, 1, 1 }, 2, 1);
		assertEquals("x^2 * y + x * y + y + x + 1.000", poly.toString());
		
		poly = new Polynomial2d(new double[] { -1, 1, 1, 1, 1, 1 }, 2, 1);
		assertEquals("-x^2 * y + x * y + y + x^2 + x + 1.000", poly.toString());
		
		poly = new Polynomial2d(new double[] { -1, 1, 1, 1, 1, -1 }, 2, 1);
		assertEquals("-x^2 * y + x * y + y + x^2 + x - 1.000", poly.toString());
		
		poly = new Polynomial2d(new double[] { -1, 1, 1, 1, 1, -1 }, 1, 2);
		assertEquals("-x * y^2 + y^2 + x * y + y + x - 1.000", poly.toString());
		
		
	}
}