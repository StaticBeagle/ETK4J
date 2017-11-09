package com.wildbitsfoundry.etk4j.util;

import java.util.Arrays;

public final class Grids {
	private Grids() {}
	
	public static class GridData {
		public final double[] X;
		public final double[] Y;
		
		public final int Rows;
		public final int Cols;
		
		private GridData(double[] x, double[] y, int rows, int cols) {
			X = x;
			Y = y;
			Rows = rows;
			Cols = cols;
		}
		
		public static GridData of(double[] x, double[] y) {
			final int xl = x.length;
			final int yl = y.length;			

			// Build grid
			x = NumArrays.repeat(x, yl);

			y = NumArrays.repeat(y, xl);
			Arrays.sort(y);
			return new GridData(x, y, yl, xl);
		}
	}
	
	public static GridData meshgrid(double[] x, double[] y) {
		return GridData.of(x, y);
	}
}
