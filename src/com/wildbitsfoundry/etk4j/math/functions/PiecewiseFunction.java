package com.wildbitsfoundry.etk4j.math.functions;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.functions.common.ExtrapolationMethod;

public abstract class PiecewiseFunction implements UnivariateFunction {
	protected double[] _x = null;
	protected double[] _coefs = null;
	ExtrapolationMethod _extrapolationMethod;

	private double _d0, _dn, _y0, _yn;

	protected PiecewiseFunction(double[] x) {
		_x = x;
		_extrapolationMethod = ExtrapolationMethod.Throw;
	}

	protected void setExtrapolationMethod(ExtrapolationMethod method) {
		_extrapolationMethod = method;
	}

	protected int findLeftIndex(double x) {
		int index = Arrays.binarySearch(_x, x);
		return index < 0 ? -(index + 2) : Math.min(index, _x.length - 2);
	}

	protected final void setEndPointsAndEndSlopes(double y0, double yn, double d0, double dn) {
		_d0 = d0;
		_dn = dn;
		_y0 = y0;
		_yn = yn;
	}

	public int getNumberOfPieces() {
		return _x.length - 1;
	}

	@Override
	public final double evaluateAt(double x) {
		if (x < _x[0]) {
			return this.extrapolate(0, x);
		}
		if (x > _x[_x.length - 1]) {
			return this.extrapolate(_x.length - 2, x);
		}
		int index = this.findLeftIndex(x);
		return this.getValueAt(index, x);
	}

	protected abstract double getValueAt(int i, double x);

	protected double extrapolate(int i, double x) {
		switch (_extrapolationMethod) {
		case ClampToZero:
			return 0.0;
		case ClampToNaN:
			return Double.NaN;
		case ClampToEndPoint:
			return i == 0 ? _y0 : _yn;
		case Natural:
			return this.getValueAt(i, x);
		case Linear:
			int n = _x.length;
			double xi = _x[0], yi = _y0, di = _d0;
			if (i != 0) {
				xi = _x[n - 1];
				yi = _yn;
				di = _dn;
			}
			return yi + di * (x - xi);
		case Periodic:
			n = _x.length;
			double delta = x - _x[0];
			double range = _x[n - 1] - _x[0];
			double xp = _x[0] + delta - range * Math.floor(delta / range);
			return this.evaluateAt(xp);
		case Throw:
			// x is smaller than the smaller value in the sequence
			if (i == 0) {
				throw new ArrayIndexOutOfBoundsException(
						String.format("x = %.4f is smaller than every number in x[]", x));
			}
			throw new ArrayIndexOutOfBoundsException(String.format("x = %.4f is bigger than every number in x[]", x));
		default:
			throw new IllegalArgumentException("invalid extrapolation option");
		}
	}
	
	protected static void checkXYDimensions(double[] x, double[] y) {
		if (x.length != y.length) {
			throw new IllegalArgumentException("x and y dimensions must be the same");
		}
	}

	protected static void checkMinkXLength(double[] x, int size) {
		if (x.length < size) {
			throw new IllegalArgumentException(String.format("x length must be >= %d", size));
		}
	}

	public double[] getBreaks() {
		return Arrays.copyOf(_x, _x.length);
	}
}
