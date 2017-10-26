package com.wildbitsfoundry.etk4j.math.interpolation2d;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.interpolation.CubicSpline;
import com.wildbitsfoundry.etk4j.math.interpolation.Interpolation;
import com.wildbitsfoundry.etk4j.math.interpolation.LinearSpline;
import com.wildbitsfoundry.etk4j.math.interpolation.Spline;

public class Spline2d {

	private double[] _y;
	private Spline[] _splines;

	private int _order;

	protected Spline2d(double[] x, double[] y, Spline[] splines, int order) {
		_y = y;
		_order = order;
		_splines = splines;
	}

	public static Spline2d newBicubicSpline(double[] x, double[] y, double[][] z) {
		final int rows = y.length;
		final int cols = x.length;
		final int order = 4;
		if (rows != cols) {
			throw new IllegalArgumentException("x and y dimensions must match");
		}
		if (rows % order != 0) {
			throw new IllegalArgumentException(String.format("y length has to be a multiple of %d", order));
		}
		if (z.length != y.length) {
			throw new IllegalArgumentException(
					String.format("The number of arrays in z has to be a multiple of %d", order));
		}
		double[] xt = Arrays.copyOf(x, cols);
		double[] yt = Arrays.copyOf(y, rows);
		Spline[] splines = new Spline[rows];
		for (int i = 0; i < rows; ++i) {
			if (cols != z[i].length) {
				throw new IllegalArgumentException("all arrays in z must have the same length");
			}
			splines[i] = CubicSpline.newCubicSplineInPlace(xt, z[i]);
		}
		return new Spline2d(xt, yt, splines, order);
	}

	public static Spline2d newBilinearSpline(double[] x, double[] y, double[][] z) {
		final int rows = y.length;
		final int cols = x.length;
		final int order = 2;
		if (rows != cols) {
			throw new IllegalArgumentException("x and y dimensions must match");
		}
		if (rows % order != 0) {
			throw new IllegalArgumentException(String.format("y length has to be a multiple of %d", order));
		}
		if (z.length != y.length) {
			throw new IllegalArgumentException(
					String.format("The number of arrays in z has to be a multiple of %d", order));
		}
		double[] xt = Arrays.copyOf(x, cols);
		double[] yt = Arrays.copyOf(y, rows);
		Spline[] splines = new Spline[rows];
		for (int i = 0; i < rows; ++i) {
			if (cols != z[i].length) {
				throw new IllegalArgumentException("all arrays in z must have the same length");
			}
			splines[i] = LinearSpline.newLinearSplineInPlace(xt, z[i]);
		}
		return new Spline2d(xt, yt, splines, order);
	}

	public double evaluateAt(double x, double y) {
		int index = this.findLeftIndex(y);

		double[] tmp = new double[_order];
		for (int i = 0; i < _order; ++i) {
			tmp[i] = _splines[i + index].evaluateAt(x);
		}
		return Interpolation.spline(Arrays.copyOfRange(_y, index, index + _order), tmp, y);
	}

	protected int findLeftIndex(double y) {
		int index = Arrays.binarySearch(_y, y);
		if(index >= 0) {
			return _order * (index / _order);
		}
		index = -(index + 2);
		boolean edge = (index + 1) % _order == 0;
		return edge ? (index + 1) - (_order >> 1) : _order * (index / _order);
	}

	public static void main(String[] args) {
		double[] x = new double[] { 1, 2, 3, 4, 5, 6, 7, 8 };
		double[] y = new double[] { 1, 2, 3, 4, 5, 6, 7, 8 };

		double[][] z = new double[][] { { 1, 4, 9, 16, 25, 36, 49, 64 }, { 4, 16, 36, 64, 100, 144, 196, 256 },
				{ 9, 36, 81, 144, 225, 324, 441, 576 }, { 16, 64, 144, 256, 400, 576, 784, 1024 },
				{ 25, 100, 225, 400, 625, 900, 1225, 1600 }, { 36, 144, 324, 576, 900, 1296, 1764, 2304 },
				{ 49, 196, 441, 784, 1225, 1764, 2401, 3136 }, { 64, 256, 576, 1024, 1600, 2304, 3136, 4096 } };

		Spline2d sp = newBilinearSpline(x, y, z);

		System.out.println(sp.evaluateAt(1, 1));
		System.out.println(sp.evaluateAt(2, 1));
		System.out.println(sp.evaluateAt(3, 1));
		System.out.println(sp.evaluateAt(4, 1));
		System.out.println(sp.evaluateAt(5, 1));
		System.out.println(sp.evaluateAt(6, 1));
		System.out.println(sp.evaluateAt(7, 1));
		System.out.println(sp.evaluateAt(8, 1));

		System.out.println(sp.evaluateAt(1, 2));
		System.out.println(sp.evaluateAt(2, 2));
		System.out.println(sp.evaluateAt(3, 2));
		System.out.println(sp.evaluateAt(4, 2));
		System.out.println(sp.evaluateAt(5, 2));
		System.out.println(sp.evaluateAt(6, 2));
		System.out.println(sp.evaluateAt(7, 2));
		System.out.println(sp.evaluateAt(8, 2));

		System.out.println(sp.evaluateAt(1, 3));
		System.out.println(sp.evaluateAt(2, 3));
		System.out.println(sp.evaluateAt(3, 3));
		System.out.println(sp.evaluateAt(4, 3));
		System.out.println(sp.evaluateAt(5, 3));
		System.out.println(sp.evaluateAt(6, 3));
		System.out.println(sp.evaluateAt(7, 3));
		System.out.println(sp.evaluateAt(8, 3));

		System.out.println(sp.evaluateAt(1, 4));
		System.out.println(sp.evaluateAt(2, 4));
		System.out.println(sp.evaluateAt(3, 4));
		System.out.println(sp.evaluateAt(4, 4));
		System.out.println(sp.evaluateAt(5, 4));
		System.out.println(sp.evaluateAt(6, 4));
		System.out.println(sp.evaluateAt(7, 4));
		System.out.println(sp.evaluateAt(8, 4));

		System.out.println(sp.evaluateAt(1, 5));
		System.out.println(sp.evaluateAt(2, 5));
		System.out.println(sp.evaluateAt(3, 5));
		System.out.println(sp.evaluateAt(4, 5));
		System.out.println(sp.evaluateAt(5, 5));
		System.out.println(sp.evaluateAt(6, 5));
		System.out.println(sp.evaluateAt(7, 5));
		System.out.println(sp.evaluateAt(8, 5));

		System.out.println(sp.evaluateAt(1, 6));
		System.out.println(sp.evaluateAt(2, 6));
		System.out.println(sp.evaluateAt(3, 6));
		System.out.println(sp.evaluateAt(4, 6));
		System.out.println(sp.evaluateAt(5, 6));
		System.out.println(sp.evaluateAt(6, 6));
		System.out.println(sp.evaluateAt(7, 6));
		System.out.println(sp.evaluateAt(8, 6));

		System.out.println(sp.evaluateAt(1, 7));
		System.out.println(sp.evaluateAt(2, 7));
		System.out.println(sp.evaluateAt(3, 7));
		System.out.println(sp.evaluateAt(4, 7));
		System.out.println(sp.evaluateAt(5, 7));
		System.out.println(sp.evaluateAt(6, 7));
		System.out.println(sp.evaluateAt(7, 7));
		System.out.println(sp.evaluateAt(8, 7));

		System.out.println(sp.evaluateAt(1, 8));
		System.out.println(sp.evaluateAt(2, 8));
		System.out.println(sp.evaluateAt(3, 8));
		System.out.println(sp.evaluateAt(4, 8));
		System.out.println(sp.evaluateAt(5, 8));
		System.out.println(sp.evaluateAt(6, 8));
		System.out.println(sp.evaluateAt(7, 8));
		System.out.println(sp.evaluateAt(8, 8));

		System.out.println(sp.evaluateAt(1.5, 1.5));

		System.out.println(sp.evaluateAt(2.5, 2.5));

		System.out.println(sp.evaluateAt(3.5, 3.5));
		
		System.out.println(sp.evaluateAt(4.5, 4.5));

		System.out.println(sp.evaluateAt(5.5, 5.5));
		
		System.out.println(sp.evaluateAt(7.5, 7.5));
	}
}
