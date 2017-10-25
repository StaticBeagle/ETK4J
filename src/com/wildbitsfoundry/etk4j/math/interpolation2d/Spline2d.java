package com.wildbitsfoundry.etk4j.math.interpolation2d;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.interpolation.CubicSpline;
import com.wildbitsfoundry.etk4j.math.interpolation.Interpolation;
import com.wildbitsfoundry.etk4j.math.interpolation.Spline;

public class Spline2d {

	private double[] _y;
	private Spline[] _splines;

	private int _numIntKnot;

	protected Spline2d(double[] x, double[] y, Spline[] splines, int knotCount) {
		_y = y;
		_numIntKnot = knotCount;
		_splines = splines;
	}

	public static Spline2d newBicubicSpline(double[] x, double[] y, double[][] z) {
		final int rows = y.length;
		final int cols = x.length;
		if (rows != cols) {
			throw new IllegalArgumentException("x and y dimensions must match");
		}
		if (z.length != y.length) {
			throw new IllegalArgumentException("not enough rows in matrix z");
		}
		double[] xt = Arrays.copyOf(x, cols);
		double[] yt = Arrays.copyOf(y, rows);
		Spline[] splines = new Spline[rows];
		for (int i = 0; i < rows; ++i) {
			if (cols != z[i].length) {
				throw new IllegalArgumentException("matrix z must be squared");
			}
			splines[i] = CubicSpline.newCubicSplineInPlace(xt, z[i]);
		}
		return new Spline2d(xt, yt, splines, 4);
	}

	public double evaluateAt(double x, double y) {
		int index = this.findLeftIndex(y);
		while (index % _numIntKnot != 0) {
			--index;
		}
		double[] tmp = new double[_numIntKnot];
		for (int i = 0; i < _numIntKnot; ++i) {
			tmp[i] = _splines[i + index].evaluateAt(x);
		}
		return Interpolation.spline(Arrays.copyOfRange(_y, index, index + _numIntKnot), tmp, y);
	}

	protected int findLeftIndex(double y) {
		int index = Arrays.binarySearch(_y, y);
		return index < 0 ? -(index + 2) : Math.min(index, _y.length - 2);
	}

	public static void main(String[] args) {
		double[] x = new double[] { 1, 2, 3, 4, 5 };
		double[] y = new double[] { 1, 2, 3, 4, 5 };

		double[][] z = new double[][] { { 1, 4, 9, 16, 25 }, { 4, 16, 36, 64, 100 }, { 9, 36, 81, 144, 225 },
				{ 16, 64, 144, 256, 400 }, { 25, 100, 225, 400, 625 } };

		Spline2d sp = newBicubicSpline(x, y, z);

		System.out.println(sp.evaluateAt(1, 1));
		System.out.println(sp.evaluateAt(2, 1));
		System.out.println(sp.evaluateAt(3, 1));
		System.out.println(sp.evaluateAt(4, 1));
		System.out.println(sp.evaluateAt(5, 1));

		System.out.println(sp.evaluateAt(1, 2));
		System.out.println(sp.evaluateAt(2, 2));
		System.out.println(sp.evaluateAt(3, 2));
		System.out.println(sp.evaluateAt(4, 2));
		System.out.println(sp.evaluateAt(5, 2));

		System.out.println(sp.evaluateAt(1, 3));
		System.out.println(sp.evaluateAt(2, 3));
		System.out.println(sp.evaluateAt(3, 3));
		System.out.println(sp.evaluateAt(4, 3));
		System.out.println(sp.evaluateAt(5, 3));

		System.out.println(sp.evaluateAt(1, 4));
		System.out.println(sp.evaluateAt(2, 4));
		System.out.println(sp.evaluateAt(3, 4));
		System.out.println(sp.evaluateAt(4, 4));
		System.out.println(sp.evaluateAt(5, 4));

		System.out.println(sp.evaluateAt(1, 5));
		System.out.println(sp.evaluateAt(2, 5));
		System.out.println(sp.evaluateAt(3, 5));
		System.out.println(sp.evaluateAt(4, 5));
		System.out.println(sp.evaluateAt(5, 5));

		System.out.println(sp.evaluateAt(1.5, 1.5));

		System.out.println(sp.evaluateAt(2.5, 2.5));

	}
}
