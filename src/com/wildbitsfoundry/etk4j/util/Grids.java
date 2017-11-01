package com.wildbitsfoundry.etk4j.util;

import java.util.Arrays;

public final class Grids {
	private Grids() {}
	
	public double[][] meshgrid(double[] x, double[] y) {
		final int xl = x.length;
		final int yl = y.length;
		
		double[][] XY = new double[2][];
		// Build grid
		x = NumArrays.repeat(x, yl);

		y = NumArrays.repeat(y, xl);
		Arrays.sort(y);
		XY[0] = x;
		XY[1] = y;
		return XY;
	}
}
