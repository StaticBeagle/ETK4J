package com.wildbitsfoundry.etk4j.math.interpolation;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.extrapolation.Extrapolators;
import com.wildbitsfoundry.etk4j.math.functions.PiecewiseFunction;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import com.wildbitsfoundry.etk4j.util.NumArrays;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkMinXLength;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkXYDimensions;

public class NearestNeighbor extends PiecewiseFunction {

	interface NeighborCalculator {
		int calculateNeighborIndex(double t);
	}
	
	private double[] _y = null;
	private double threshold = 0.5;
	private NeighborCalculator neighborCalculator;

	protected NearestNeighbor(double[] x, double[] y, double threshold, boolean roundUp) {
		super(x);
		_y = Arrays.copyOf(y, y.length);
		this.threshold = threshold;
		if(roundUp) {
			this.neighborCalculator = new NeighborCalculator() {
				@Override
				public int calculateNeighborIndex(double t) {
					return t >= threshold ? 1 : 0;
				}
			};
		} else {
			this.neighborCalculator = new NeighborCalculator() {
				@Override
				public int calculateNeighborIndex(double t) {
					return t > threshold ? 1 : 0;
				}
			};
		}
		if(!NumArrays.isAscending(x)) {
			throw new IllegalArgumentException("x must be monotonically increasing");
		}
		double x0 = this.x[0];
		double xn = this.x[this.x.length - 1];
		double y0 = this.evaluateAt(x0);
		double yn = this.evaluateAt(xn);
		this.setExtrapolator(new Extrapolators.ClampToEndPointExtrapolator(x0, xn, y0, yn));
	}
	
	public static NearestNeighbor newNearestNeighbor(double[] x, double[] y) {
		return newNearestNeighborInPlace(Arrays.copyOf(x, x.length), y, 0.5, true);
	}

	public static NearestNeighbor newNearestNeighbor(double[] x, double[] y, double threshold) {
		return newNearestNeighborInPlace(Arrays.copyOf(x, x.length), y, threshold, true);
	}

	public static NearestNeighbor newNearestNeighbor(double[] x, double[] y, boolean roundUp) {
		return newNearestNeighborInPlace(Arrays.copyOf(x, x.length), y, 0.5, roundUp);
	}

	public static NearestNeighbor newNearestNeighbor(double[] x, double[] y, double threshold, boolean roundUp) {
		return newNearestNeighborInPlace(Arrays.copyOf(x, x.length), y, threshold, roundUp);
	}
	
	public static NearestNeighbor newNearestNeighborInPlace(double[] x, double[] y) {
		return newNearestNeighborInPlace(x, y, 0.5, true);
	}

	public static NearestNeighbor newNearestNeighborInPlace(double[] x, double[] y, double threshold) {
		return newNearestNeighborInPlace(x, y, threshold, true);
	}

	public static NearestNeighbor newNearestNeighborInPlace(double[] x, double[] y, boolean roundUp) {
		return newNearestNeighborInPlace(x, y, 0.5, roundUp);
	}

	public static NearestNeighbor newNearestNeighborInPlace(double[] x, double[] y, double threshold, boolean roundUp) {
		checkXYDimensions(x, y);
		checkMinXLength(x, 2);
		if(threshold < 0 || threshold > 1) {
			throw new IllegalArgumentException("Threshold must be 0 < threshold <= 1.");
		}

		NearestNeighbor nh = new NearestNeighbor(x, y, threshold, roundUp);
		return nh;
	}


	@Override
	public double evaluateAt(int index, double x) {
		double t = (x - this.x[index]) / (this.x[index + 1] - this.x[index]);
		return _y[index + neighborCalculator.calculateNeighborIndex(t)];
	}

	@Override
	public UnivariateFunction getSegment(int index) {
		final double yi = this.evaluateAt(index, x[index]);
		UnivariateFunction fn = new UnivariateFunction() {
			
			@Override
			public double evaluateAt(double x) {
				return yi;
			}
		};
		return fn;
	}

	public static void main(String[] args) {
		double[] x = { 1, 2, 3, 4 };
		double[] y = { 1, 4, 9, 16 };

		NearestNeighbor nh = NearestNeighbor.newNearestNeighbor(x, y);
		double[] xi = {1, 1.5, 2, 2.5, 3, 3.5, 4};
		System.out.println();
		for(double d : xi) {
			System.out.printf("x: %f, y: %f%n", d, nh.evaluateAt(d));
		}
	}
}
