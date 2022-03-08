package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.math.extrapolation.ExtrapolationMethod;
import com.wildbitsfoundry.etk4j.math.extrapolation.Extrapolators;
import com.wildbitsfoundry.etk4j.math.functions.DifferentiableFunction;
import com.wildbitsfoundry.etk4j.math.functions.IntegrableFunction;
import com.wildbitsfoundry.etk4j.math.functions.PiecewiseFunction;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import com.wildbitsfoundry.etk4j.util.DoubleArrays;

import java.util.Arrays;


public abstract class Spline extends PiecewiseFunction implements DifferentiableFunction, IntegrableFunction {
	
	private int order;
	private double[] _indefiniteIntegral = null;

	protected Spline(double[] x, int order) {
		super(x);
		this.order = order;
	}
	
	@Override
	public UnivariateFunction getSegment(int index) {
		int offset = index * order;
		final double[] coefficients = Arrays.copyOfRange(this.coefficients, offset, offset + order);
		final double x0 = x[index];
		UnivariateFunction fn = x -> DoubleArrays.horner(coefficients, x - x0);
		return fn;
	}

	protected double evaluateDerivativeAt(int index, double x) {
		double t = x - this.x[index];
		index *= order;
		final int length = order - 1;
		double result = 0.0;
		for (int j = 0; j < length; ++j) {
			result *= t;
			result += coefficients[index++] * (length - j);
		}
		return result;
	}

	public double evaluateAntiDerivativeAt(int index, double x) {
		double t = x - this.x[index];
		index *= order;
		double result = 0.0;
		for (int j = 0; j < order; ++j) {
			result += coefficients[index++] / (order - j);
			result *= t;
		}
		return result;
	}
	
	@Override
	public double evaluateAt(int index, double x) {
		double t = x - this.x[index];
		index *= order;
		double result = 0;
		for (int j = 0; j < order; ++j) {
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
		if(b < a) {
			throw new IllegalArgumentException("The upper integration limit b has to be greater or equal to the lower" +
					" integration limit a.");
		}
		if (a < x[0] || b > x[x.length - 1]) {
			throw new IllegalArgumentException(
					String.format("The spline is not defined outside of [%.4g, %.4g]", x[0], x[x.length - 1]));
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

		final int size = coefficients.length / order;
		_indefiniteIntegral = new double[size];
		for (int i = 0; i < size - 1; ++i) {
			_indefiniteIntegral[i + 1] = _indefiniteIntegral[i] + this.evaluateAntiDerivativeAt(i, x[i + 1]);
		}
	}

	public void setExtrapolationMethod(ExtrapolationMethod method) {
		double x0 = x[0];
		double xn = x[x.length - 1];
		double y0 = this.evaluateAt(x0);
		double yn = this.evaluateAt(xn);
		switch (method) {
			case CLAMP_TO_ZERO:
				this.setExtrapolator(new Extrapolators.ClampToZeroExtrapolator());
				break;
			case CLAMP_TO_NAN:
				this.setExtrapolator(new Extrapolators.ClampToNaNExtrapolator());
				break;
			case CLAMP_TO_END_POINT:
				this.setExtrapolator(new Extrapolators.ClampToEndPointExtrapolator(x0, xn, y0, yn));
				break;
			case NATURAL:
				UnivariateFunction lfn = this.getFirstSegment();
				UnivariateFunction rfn = this.getLastSegment();
				this.setExtrapolator(new Extrapolators.NaturalExtrapolator(lfn, rfn, x0, xn));
				break;
			case LINEAR:
				this.setExtrapolator(new Extrapolators.LinearExtrapolator(this, x0, xn, y0, yn));
				break;
			case THROW:
				this.setExtrapolator(new Extrapolators.ThrowExtrapolator(x0, xn));
				break;
			default:
				throw new IllegalArgumentException("Invalid extrapolation option.");
		}
	}
}
