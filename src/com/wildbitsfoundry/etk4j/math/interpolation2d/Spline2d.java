package com.wildbitsfoundry.etk4j.math.interpolation2d;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.interpolation.CubicSpline;
import com.wildbitsfoundry.etk4j.math.interpolation.Interpolation;
import com.wildbitsfoundry.etk4j.math.interpolation.Spline;

public class Spline2d {
	
	private double[] _x;
	private double[] _y;
	private Spline[] _splines;
	
	protected Spline2d(double[] x, double[] y, double[][] z) {
		// check x and y dims match
		if(z.length != y.length) {
			throw new IllegalArgumentException("not enough rows in matrix z");
		}
		
		// assuming type cubic
		for(int i = 0; i < z.length; ++i) {
			double[] dydx = CubicSpline.buildNotAKnotSpline(x, z[i]);
			_splines[i] = CubicSpline.newHermiteSplineInPlace(x, y, dydx);
		}
		
	}
	
	public double evaluateAt(double x, double y) {
		int index = this.findLeftIndex(y);
		int order = 4;
		while(index % order != 0) {
			--index;
		}
		double[] tmp = new double[order];
		for(int i = 0; i < order; ++i) {
			tmp[i] = _splines[i + index].evaluateAt(x);
		}
		
		double[] gg = Arrays.copyOfRange(_y, index, index + order);
		return Interpolation.spline(gg, tmp, y);
		
		
	}
	
	protected int findLeftIndex(double y) {
		int index = Arrays.binarySearch(_y, y);
		return index < 0 ? -(index + 2) : Math.min(index, _y.length - 2);
	}
}
