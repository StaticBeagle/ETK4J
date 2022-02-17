package com.wildbitsfoundry.etk4j.math.interpolation;

import java.util.Arrays;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkMinXLength;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkXYDimensions;

public class LinearSpline extends Spline {

	protected LinearSpline(double[] x, double[] y) {
		super(x, 2);

		final int n = this.x.length;
		// compute coefficients
		coefficients = new double[(n - 1) * 2]; // 2 coefficients and n - 1 segments
		for (int i = 0, j = 0; i < n - 1; ++i, ++j) {
			double hx = this.x[i + 1] - this.x[i];
			double a = (y[i + 1] - y[i]) / hx;
			double b = y[i];
			coefficients[j] = a;
			coefficients[++j] = b;
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
	public double evaluateAt(int i, double x) {
		double t = x - this.x[i];
		i <<= 1;
		return coefficients[i + 1] + t * coefficients[i];
	}

	@Override
	protected double evaluateDerivativeAt(int i, double x) {
		i <<= 1;
		return coefficients[i];
	}

	@Override
	public double evaluateAntiDerivativeAt(int i, double x) {
		double t = x - this.x[i];
		i <<= 1;
		return t * (coefficients[i + 1] + t * coefficients[i] * 0.5);
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0, j = 0; i < x.length - 1; ++i, ++j) {
			double a = coefficients[j];
			double b = coefficients[++j];

			sb.append("S").append(i + 1).append("(x) = ")
					.append(a != 0d ? String.format("%.4g * (x - %.4f)", a, x[i]) : "")
					.append(b != 0d ? String.format(" + %.4g", b) : "")
					.append(System.lineSeparator());
		}
		sb.setLength(Math.max(sb.length() - System.lineSeparator().length(), 0));
		return sb.toString().replace("+ -", "- ").replace("- -", "+ ")
				.replace("=  + ", "= ").replace("=  - ", "= -");
	}
}
