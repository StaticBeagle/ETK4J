package com.wildbitsfoundry.etk4j.math.functions;

import java.util.Arrays;

public abstract class PiecewiseFunction implements UnivariateFunction {
	
	public enum ExtrapolationMethod {
		ClampToZero, ClampToNaN, ClampToEndPoint, Natural, Linear, Periodic, Throw
	}

	protected double[] _x = null;
	protected final double _y0, _yn;
	private double _d0 = 0.0, _dn = 0.0;
	private ExtrapolationMethod _extrapolationMethod = ExtrapolationMethod.Throw;
	
	protected PiecewiseFunction(double[] x, double y0, double yn) {
		_x = x;
		_y0 = y0;
		_yn = yn;
	}
	
	protected void setExtrapolationMethod(ExtrapolationMethod method) {
		_extrapolationMethod = method;
	}

	protected int findIndex(double x) {
		int index = Arrays.binarySearch(_x, x);
		return index < 0 ? -(index + 2) : Math.min(index, _x.length - 2);
	}
	
	protected final void setEndSlopes(double d0, double dn) {
		_d0 = d0;
		_dn = dn;
	}
	
	public int getSize() {
		return _x.length - 1;
	}

	@Override
	public final double evaluateAt(double x) {
		if(x < _x[0]) {
			return this.extrapolate(0, x);
		}
		if(x > _x[_x.length - 1]) {
			return this.extrapolate(_x.length - 2, x);
		}
		int index = this.findIndex(x);
		return this.getValueAt(index, x);
	}
	
	//public abstract <T extends UnivariateFunction> T getSegment(int index); 

	protected abstract double getValueAt(int index, double x);
	
	protected double extrapolate(int index, double x) {
		switch (_extrapolationMethod) {
		case ClampToZero:
			return 0.0;
		case ClampToNaN:
			return Double.NaN;
		case ClampToEndPoint:
			return index == 0 ? _y0 : _yn;
		case Natural:
			return this.getValueAt(index, x);
		case Linear:
			int n = _x.length;
			double xi = _x[0], yi = _y0, di = _d0;
			if (index != 0) {
				xi = _x[n - 1];
				yi = _yn;
				di =_dn;
			}
			return yi + di * (x - xi);
		case Periodic:
			n = _x.length;
			double delta = x - _x[0];
			double range = _x[n - 1] - _x[0];
			double xp = _x[0] + delta - range * Math.floor(delta / range);
			return this.evaluateAt(xp);
		default:
			// x is smaller than the smaller value in the sequence
			if (index == 0) {
				throw new ArrayIndexOutOfBoundsException(
						String.format("x = %.4f is smaller than every number in x[]", x));
			}
			throw new ArrayIndexOutOfBoundsException(
					String.format("x = %.4f is bigger than every number in x[]", x));
		}
	}
}
