package com.wildbitsfoundry.etk4j.util;

import java.util.Arrays;

public final class Grids {
	private Grids() {}
	
	public static class GridData {
		public double[] X;
		public double[] Y;
		
		private int _rows;
		private int _cols;
		
		private GridData() {}
		
		public static GridData of(double[] x, double[] y) {
			final int xl = x.length;
			final int yl = y.length;
			
			GridData gd = new GridData();
			// Build grid
			x = NumArrays.repeat(x, yl);

			y = NumArrays.repeat(y, xl);
			Arrays.sort(y);
			gd._rows = yl;
			gd._cols = xl;
			gd.X = x;
			gd.Y = y;
			return gd;
		}
		
		public int getRowCount() {
			return _rows;
		}
		
		public int getColumnCount() {
			return _cols;
		}
	}
	
	public static GridData meshgrid2d(double[] x, double[] y) {
		return GridData.of(x, y);
	}
}
