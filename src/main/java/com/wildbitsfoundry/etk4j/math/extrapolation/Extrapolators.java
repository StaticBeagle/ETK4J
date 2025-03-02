package com.wildbitsfoundry.etk4j.math.extrapolation;

import com.wildbitsfoundry.etk4j.math.function.DifferentiableFunction;
import com.wildbitsfoundry.etk4j.math.function.UnivariateFunction;

public final class Extrapolators {

	private Extrapolators() {
	}

	/**
	 * The {@code ClampToEndPointExtrapolator} class clamps the value of the extrapolation to the values prescribed by
	 * {@code y0 and yn}.
	 */
	public static class ClampToEndPointExtrapolator implements Extrapolator {

		private double x0, xn, y0, yn;

		public ClampToEndPointExtrapolator(double x0, double xn, double y0, double yn) {
			this.x0 = x0;
			this.xn = xn;
			this.y0 = y0;
			this.yn = yn;
		}

		/**
		 * Extrapolation.
		 * @param x Argument at which to evaluate the extrapolation.
		 * @return y0 if x &lt; x0 and yn if x &gt; xn.
		 */
		@Override
		public double extrapolate(double x) {
			if (x < x0) {
				return y0;
			} else if (x > xn) {
				return yn;
			}
			throw new IllegalArgumentException(String.format("%.4f is not outside of [%.4f, %.4f]", x, x0, xn));
		}
	}

	/**
	 * The {@code ClampToValueExtrapolator} class clamps the value of the extrapolation to a prescribed value.
	 */
	public static class ClampToValueExtrapolator implements Extrapolator {

		private double value;

		public ClampToValueExtrapolator(double value) {
			this.value = value;
		}

		/**
		 * Extrapolate.
		 * @param x Argument at which to evaluate the extrapolation.
		 * @return The value assigned during construction.
		 */
		@Override
		public double extrapolate(double x) {
			return value;
		}

	}

	/**
	 * The {@code ClampToNaNExtrapolator} class returns NaN as the result of the extrapolation.
	 */
	public static class ClampToNaNExtrapolator extends ClampToValueExtrapolator {

		public ClampToNaNExtrapolator() {
			super(Double.NaN);
		}
	}

	/**
	 * The {@code ClampToNaNExtrapolator} class returns zero as the result of the extrapolation.
	 */
	public static class ClampToZeroExtrapolator extends ClampToValueExtrapolator {

		public ClampToZeroExtrapolator() {
			super(0.0);
		}
	}

	/**
	 * The {@code LinearExtrapolator} class constructs a line {@code y = m * (x - xi) + yi} where {@code m} is calculated by
	 * differentiating the {@link DifferentiableFunction}. {@code (xi, yi)} are the left-hand side coordinates if x is
	 * less than the minimum value for extrapolation or the right-hand side values otherwise.
	 */
	public static class LinearExtrapolator implements Extrapolator {

		private DifferentiableFunction dfn;
		private double x0, xn, y0, yn;

		public LinearExtrapolator(DifferentiableFunction dfn, double x0, double xn, double y0, double yn) {
			this.dfn = dfn;
			this.x0 = x0;
			this.xn = xn;
			this.y0 = y0;
			this.yn = yn;
		}

		@Override
		public double extrapolate(double x) {
			double xi = x0;
			double yi = y0;
			if (x > xn) {
				xi = xn;
				yi = yn;
			}
			double di = dfn.differentiate(xi);
			return yi + di * (x - xi);
		}
	}

	/**
	 * The {@code NaturalExtrapolator} class evaluates the extrapolation at a function specified for the left-hand side
	 * of the extrapolation if x is on the left-hand side of the extrapolation or evaluates the extrapolation at a
	 * function specified for the right-hand side otherwise.
	 */
	public static class NaturalExtrapolator implements Extrapolator {

		private UnivariateFunction leftFn;
		private UnivariateFunction rightFn;

		private double x0, xn;

		public NaturalExtrapolator(UnivariateFunction leftFn, UnivariateFunction rightFn, double x0, double xn) {
			this.leftFn = leftFn;
			this.rightFn = rightFn;
			this.x0 = x0;
			this.xn = xn;
		}

		@Override
		public double extrapolate(double x) {
			if (x < x0) {
				return leftFn.evaluateAt(x);
			} else if (x > xn) {
				return rightFn.evaluateAt(x);
			}
			throw new IllegalArgumentException(String.format("%.4f is not outside of [%.4f, %.4f]", x, x0, xn));
		}
	}

	/**
	 * The {@code ThrowExtrapolator} class throws and {@link ArrayIndexOutOfBoundsException} if x is outside of the
	 * valid range or {@link IllegalArgumentException} if x is within the allowed range.
	 */
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
