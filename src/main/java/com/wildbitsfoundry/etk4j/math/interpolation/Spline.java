package com.wildbitsfoundry.etk4j.math.interpolation;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.functions.DifferentiableFunction;
import com.wildbitsfoundry.etk4j.math.functions.IntegrableFunction;
import com.wildbitsfoundry.etk4j.math.functions.PiecewiseFunction;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import com.wildbitsfoundry.etk4j.math.interpolation.Extrapolators.*;
import com.wildbitsfoundry.etk4j.util.NumArrays;


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
		final double[] coefs = Arrays.copyOfRange(coefficients, offset, offset + _order);
		final double x0 = _x[index];
		UnivariateFunction fn = x -> NumArrays.horner(coefs, x - x0);
		return fn;
	}

	// TODO shouldn't this be in the parent?
	public void setExtrapolationMethod(ExtrapolationMethod method) {
		double x0 = _x[0];
		double xn = _x[_x.length - 1];
		double y0 = this.evaluateAt(x0);
		double yn = this.evaluateAt(xn);
		switch (method) {
		case CLAMP_TO_ZERO:
			this.setExtrapolator(new ClampToZeroExtrapolator());
			break;
		case CLAMP_TO_NAN:
			this.setExtrapolator(new ClampToNaNExtrapolator());
			break;
		case CLAMP_TO_END_POINT:
			this.setExtrapolator(new ClampToEndPointExtrapolator(x0, xn, y0, yn));
			break;
		case NATURAL:
			UnivariateFunction lfn = this.getFirstSegment();
			UnivariateFunction rfn = this.getLastSegment();
			this.setExtrapolator(new NaturalExtrapolator(lfn, rfn, x0, xn));
			break;
		case LINEAR:
			this.setExtrapolator(new LinearExtrapolator(this, x0, xn, y0, yn));
			break;
		case PERIODIC:
			this.setExtrapolator(new PeriodicExtrapolator(this, x0, xn));
			break;
		case THROW:
			this.setExtrapolator(new ThrowExtrapolator(x0, xn));
			break;
		default:
			throw new IllegalArgumentException("Invalid extrapolation option");
		}
	}

	// TODO this t and x is confusing, the interface should expose x not t
	protected double evaluateDerivativeAt(int index, double x) {
		double t = x - _x[index];
		index *= _order;
		final int length = _order - 1;
		double result = 0.0;
		for (int j = 0; j < length; ++j) {
			result *= t;
			result += coefficients[index++] * (length - j);
		}
		return result;
	}

	// TODO this t and x is confusing, the interface should expose x not t
	public double evaluateAntiDerivativeAt(int index, double x) {
		double t = x - _x[index];
		index *= _order;
		double result = 0.0;
		for (int j = 0; j < _order; ++j) {
			result += coefficients[index++] / (_order - j);
			result *= t;
		}
		return result;
	}
	
	@Override
	public double evaluateAt(int index, double x) {
		double t = x - _x[index];
		index *= _order;
		double result = 0;
		for (int j = 0; j < _order; ++j) {
			result *= t;
			result += coefficients[index++];
		}
		return result;
	}


	@Override
	public double differentiate(double x) {
		int i = this.findSegmentIndex(x);
		return this.evaluateDerivativeAt(i, x);
	}

	@Override
	public double integrate(double a, double b) {
		if (a < _x[0] || b > _x[_x.length - 1]) {
			throw new IllegalArgumentException(
					String.format("The spline is not defined outside of [%.4g, %.4g]", _x[0], _x[_x.length - 1]));
		}
		return integrate(b) - integrate(a);
	}

	public double integrate(double x) {
		// Lazy creating of values
		if (_indefiniteIntegral == null) {
			this.calculateIntegral();
		}

		int i = this.findSegmentIndex(x);
		return _indefiniteIntegral[i] + this.evaluateAntiDerivativeAt(i, x);
	}

	private void calculateIntegral() {
		// Lazy creating the values
		if (_indefiniteIntegral != null) {
			return;
		}

		final int size = coefficients.length / _order;
		_indefiniteIntegral = new double[size];
		for (int i = 0; i < size - 1; ++i) {
			_indefiniteIntegral[i + 1] = _indefiniteIntegral[i] + this.evaluateAntiDerivativeAt(i, _x[i + 1]);
		}
	}
}
