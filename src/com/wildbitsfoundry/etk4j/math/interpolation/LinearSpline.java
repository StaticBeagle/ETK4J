package com.wildbitsfoundry.etk4j.math.interpolation;

import java.util.Arrays;

public class LinearSpline extends Spline {
	

	protected LinearSpline(double[] x, double[] y) {
		super(x, 1, y[0], y[y.length - 1]);
		
		final int n = _x.length;
		// compute coefficients
		_coefs = new double[(n - 1) * 2]; // 2 coefficients and n - 1 segments
		for(int i = 0, j = 0; i < n - 1; ++i, ++j) {
			double hx = _x[i + 1] - _x[i];
			double a = (y[i + 1] - y[i]) / hx;
			double b = y[i];
			_coefs[j] = a;
			_coefs[++j] = b;
		}
		double d0 = _coefs[0];
		double dn = _coefs[_coefs.length - 2];
		super.setEndSlopes(d0, dn);
	}
	
	public static LinearSpline newLinearSpline(double[] x, double[] y) {
		return newLinearSplineInPlace(Arrays.copyOf(x, x.length), y);
	}
	
	public static LinearSpline newLinearSplineInPlace(double[] x, double[] y) {
		checkXYDimensions(x, y);
		checkMinkXLength(x, 2);
		return new LinearSpline(Arrays.copyOf(x, x.length), y);
	}
	
	@Override
	public void setExtrapolationMethod(ExtrapolationMethod method) {
		super.setExtrapolationMethod(method);
	}

	@Override
	public double differentiate(double x) {
		throw new RuntimeException("Method not implemented yet");
	}

	
	@Override
	protected double evaluateAntiDerivativeAt(int i, double t) {
		i = i << 1;
		return t * (_coefs[i + 1] + t * _coefs[i] * 0.5);
	}
	
	public static void main(String[] args) {
		double[] x = new double[] {1, 2, 3, 4};
		double[] y = new double[] {1, 4, 9, 16};
		
		LinearSpline ls = newLinearSpline(x, y);
		CubicSpline cs = CubicSpline.newNotAKnotSpline(x, y);
		NearestNeighbor nh = NearestNeighbor.newNearestNeighbor(x, y);
		
		System.out.println(ls.evaluateAt(4));
		
		System.out.println(ls.integrate(1, 4));
		
		System.out.println(cs.evaluateAt(4));
		
		System.out.println(cs.integrate(1, 4));
		
		System.out.println(nh.evaluateAt(1.6));
	}
}
