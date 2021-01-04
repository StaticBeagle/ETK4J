package com.wildbitsfoundry.etk4j.math.interpolation;

import java.util.Arrays;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkMinXLength;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkXYDimensions;

public class LinearSpline extends Spline {

	protected LinearSpline(double[] x, double[] y) {
		super(x, 2);

		final int n = _x.length;
		// compute coefficients
		_coefs = new double[(n - 1) * 2]; // 2 coefficients and n - 1 segments
		for (int i = 0, j = 0; i < n - 1; ++i, ++j) {
			double hx = _x[i + 1] - _x[i];
			if (hx <= 0.0) {
				throw new IllegalArgumentException("x must be monotonically increasing");
			}
			double a = (y[i + 1] - y[i]) / hx;
			double b = y[i];
			_coefs[j] = a;
			_coefs[++j] = b;
		}
	}

	public static LinearSpline newLinearSpline(double[] x, double[] y) {
		return newLinearSplineInPlace(Arrays.copyOf(x, x.length), y);
	}

	public static LinearSpline newLinearSplineInPlace(double[] x, double[] y) {
		checkXYDimensions(x, y);
		checkMinXLength(x, 2);
		return new LinearSpline(x, y);
	}

	@Override
	public double evaluateSegmentAt(int i, double x) {
		double t = x - _x[i];
		i <<= 1;
		return _coefs[i + 1] + t * _coefs[i];
	}

	@Override
	protected double evaluateDerivativeAt(int i, double t) {
		i = i << 1;
		return _coefs[i];
	}

	@Override
	public double evaluateAntiDerivativeAt(int i, double t) {
		i = i << 1;
		return t * (_coefs[i + 1] + t * _coefs[i] * 0.5);
	}

	public static void main(String[] args) {
		double[] x = { 1, 2, 3, 4 };
		double[] y = { 1, 4, 9, 16 };

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
