package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.math.extrapolation.Extrapolator;
import com.wildbitsfoundry.etk4j.math.functions.DifferentiableFunction;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;

final class Extrapolators {

	private Extrapolators() {

	}

	public static class ClampToEndPointExtrapolator implements Extrapolator {

		private double _x0, _xn, _y0, _yn;

		public ClampToEndPointExtrapolator(double x0, double xn, double y0, double yn) {
			_x0 = x0;
			_xn = xn;
			_y0 = y0;
			_yn = yn;
		}

		@Override
		public double extrapolate(double x) {
			if (x < _x0) {
				return _y0;
			} else if (x > _xn) {
				return _yn;
			}
			throw new IllegalArgumentException(String.format("%.4f is not outside of [%.4f, %.4f]", x, _x0, _xn));
		}
	}

	public static class ClampToValueExtrapolator implements Extrapolator {

		private double _value;

		public ClampToValueExtrapolator(double value) {
			_value = value;
		}

		@Override
		public double extrapolate(double x) {
			return _value;
		}

	}

	public static class ClampToNaNExtrapolator extends ClampToValueExtrapolator {

		public ClampToNaNExtrapolator() {
			super(Double.NaN);
		}
	}

	public static class ClampToZeroExtrapolator extends ClampToValueExtrapolator {

		public ClampToZeroExtrapolator() {
			super(0.0);
		}
	}

	public static class LinearExtrapolator implements Extrapolator {

		private DifferentiableFunction _dfn;
		private double _x0, _xn, _y0, _yn;

		public LinearExtrapolator(DifferentiableFunction dfn, double x0, double xn, double y0, double yn) {
			_dfn = dfn;
			_x0 = x0;
			_xn = xn;
			_y0 = y0;
			_yn = yn;
		}

		@Override
		public double extrapolate(double x) {
			double xi = _x0;
			double yi = _y0;
			if (x > _xn) {
				xi = _xn;
				yi = _yn;
			}
			double di = _dfn.differentiate(xi);
			return yi + di * (x - xi);
		}
	}

	public static class NaturalExtrapolator implements Extrapolator {

		private UnivariateFunction _leftFn;
		private UnivariateFunction _rightFn;

		private double _x0, _xn;

		public NaturalExtrapolator(UnivariateFunction leftFn, UnivariateFunction rightFn, double x0, double xn) {
			_leftFn = leftFn;
			_rightFn = rightFn;
			_x0 = x0;
			_xn = xn;
		}

		@Override
		public double extrapolate(double x) {
			if (x < _x0) {
				return _leftFn.evaluateAt(x);
			} else if (x > _xn) {
				return _rightFn.evaluateAt(x);
			}
			throw new IllegalArgumentException(String.format("%.4f is not outside of [%.4f, %.4f]", x, _x0, _xn));
		}
	}
	
	public static class PeriodicExtrapolator implements Extrapolator {

		private UnivariateFunction _fn;
		private double _x0, _xn;
		
		public PeriodicExtrapolator(UnivariateFunction fn, double x0, double xn) {
			_fn = fn;
			_x0 = x0;
			_xn = xn;
		}

		@Override
		public double extrapolate(double x) {
			double delta = x - _x0;
			double range = _xn- _x0;
			double xp = _x0 + delta - range * Math.floor(delta / range);
			return _fn.evaluateAt(xp);
		}
	}
	
	public static class ThrowExtrapolator implements Extrapolator {
		
		private double _x0, _xn;
		
		public ThrowExtrapolator(double x0, double xn) {
			_x0 = x0;
			_xn = xn;
		}

		@Override
		public double extrapolate(double x) {
			// x is smaller than the smaller value in the sequence
			if (x < _x0) {
				throw new ArrayIndexOutOfBoundsException(
						String.format("x = %.4f is smaller than every number in x[]", x));
			} else if(x > _xn) {
				throw new ArrayIndexOutOfBoundsException(String.format("x = %.4f is bigger than every number in x[]", x));			
			}
			throw new IllegalArgumentException(String.format("%.4f is not outside of [%.4f, %.4f]", x, _x0, _xn));
		}
	}
}
