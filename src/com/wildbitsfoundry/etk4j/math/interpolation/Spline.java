package com.wildbitsfoundry.etk4j.math.interpolation;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.functions.DifferentiableFunction;
import com.wildbitsfoundry.etk4j.math.functions.IntegrableFunction;
import com.wildbitsfoundry.etk4j.math.functions.PiecewiseFunction;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import com.wildbitsfoundry.etk4j.math.interpolation.Extrapolators.*;


public abstract class Spline extends PiecewiseFunction implements DifferentiableFunction, IntegrableFunction {
	
	private int _order;
	private double[] _indefiniteIntegral = null;

	protected Spline(double[] x, int order) {
		super(x);
		_order = order;
	}
	
	@Override
	public UnivariateFunction getSegment(int index) {
		int offset = index * _order;
		final double[] coefs = Arrays.copyOfRange(_coefs, offset, offset + _order);
		final double x0 = _x[index];
		UnivariateFunction fn = new UnivariateFunction() {
			
			@Override
			public double evaluateAt(double x) {
				double t = x - x0;
				double result = 0.0;
				for (int j = 0; j < _order; ++j) {
					result = result * t + coefs[j];
				}
				return result;
			}
		};
		return fn;
	}
	
	public void setExtrapolationMethod(ExtrapolationMethod method) {
		double x0 = _x[0];
		double xn = _x[_x.length - 1];
		double y0 = this.evaluateAt(x0);
		double yn = this.evaluateAt(xn);
		switch (method) {
		case ClampToZero:
			this.setExtrapolator(new ClampToZeroExtrapolator());
			break;
		case ClampToNaN:
			this.setExtrapolator(new ClampToNaNExtrapolator());
			break;
		case ClampToEndPoint:
			this.setExtrapolator(new ClampToEndPointExtrapolator(x0, xn, y0, yn));
			break;
		case Natural:
			UnivariateFunction lfn = this.getFirstSegment();
			UnivariateFunction rfn = this.getLastSegment();
			this.setExtrapolator(new NaturalExtrapolator(lfn, rfn, x0, xn));
			break;
		case Linear:
			this.setExtrapolator(new LinearExtrapolator(this, x0, xn, y0, yn));
			break;
		case Periodic:
			this.setExtrapolator(new PeriodicExtrapolator(this, x0, xn));
			break;
		case Throw:
			this.setExtrapolator(new ThrowExtrapolator(x0, xn));
			break;
		default:
			throw new IllegalArgumentException("invalid extrapolation option");
		}
	}
		
	protected abstract double evaluateDerivativeAt(int i, double t);
	protected abstract double evaluateAntiDerivativeAt(int i, double t);
	
	@Override
	public double evaluateSegmentAt(int i, double x) {

		double t = x - _x[i];
		i *= _order;
		double result = 0;
		for (int j = 0; j < _order; ++j) {
			result = result * t + _coefs[i++];
		}
		return result;
	}


	@Override
	public double differentiate(double x) {
		int i = this.findSegmentIndex(x);
		double t = x - _x[i];
		return this.evaluateDerivativeAt(i, t);
	}

	@Override
	public double integrate(double a, double b) {
		if (a < _x[0] || b > _x[_x.length - 1]) {
			throw new IllegalArgumentException(
					String.format("The spline is not defined outside of [%.4g, %.4g]", _x[0], _x[_x.length - 1]));
		}
		return integrate(b) - integrate(a);
	}

	private double integrate(double x) {
		// Lazy creating of values
		if (_indefiniteIntegral == null) {
			this.calcuateIntegral();
		}

		int i = this.findSegmentIndex(x);
		double t = x - _x[i];
		return _indefiniteIntegral[i] + this.evaluateAntiDerivativeAt(i, t);
	}

	private void calcuateIntegral() {
		// Lazy creating the values
		if (_indefiniteIntegral != null) {
			return;
		}

		final int size = _coefs.length / _order;
		_indefiniteIntegral = new double[size];
		for (int i = 0; i < size - 1; ++i) {
			double t = _x[i + 1] - _x[i];
			_indefiniteIntegral[i + 1] = _indefiniteIntegral[i] + this.evaluateAntiDerivativeAt(i, t);
		}
	}
}
