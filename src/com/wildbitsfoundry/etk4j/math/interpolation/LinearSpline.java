package com.wildbitsfoundry.etk4j.math.interpolation;

public class LinearSpline extends Spline {
	
	private double[] _coefs = null;

	protected LinearSpline(double[] x, double[] y) {
		super(x, y[0], y[y.length - 1]);
		
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
		checkXYDimensions(x, y);
		checkMinkXLength(x, 2);
		return new LinearSpline(x, y);
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
	public double integrate(double x0, double x1) {
		throw new RuntimeException("Method not implemented yet");
	}

	@Override
	protected double getValueAt(int index, double x) {
		double t = x - _x[index];
		index = index << 1;
		return _coefs[index + 1] + t * _coefs[index];
	}
}
