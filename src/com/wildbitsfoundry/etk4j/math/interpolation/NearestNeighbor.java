package com.wildbitsfoundry.etk4j.math.interpolation;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.functions.PiecewiseFunction;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import com.wildbitsfoundry.etk4j.util.NumArrays;
import static com.wildbitsfoundry.etk4j.math.interpolation.Extrapolators.ClampToEndPointExtrapolator;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkMinXLength;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkXYDimensions;

public class NearestNeighbor extends PiecewiseFunction {
	
	interface NeighborCalculator {
		int calculateNeighborIndex(double t);
	}
	
	private double[] _y = null;
	private NeighborCalculator _calculator = null;
	protected NearestNeighbor(double[] x, double[] y) {
		super(x);
		_y = Arrays.copyOf(y, y.length);
		if(!NumArrays.isAscending(x)) {
			throw new IllegalArgumentException("x must be monotonically increasing");
		}
		double x0 = _x[0];
		double xn = _x[_x.length - 1];
		double y0 = this.evaluateAt(x0);
		double yn = this.evaluateAt(xn);
		this.setExtrapolator(new ClampToEndPointExtrapolator(x0, xn, y0, yn));
	}
	
	public static NearestNeighbor newNearestNeighbor(double[] x, double[] y) {
		return newNearestNeighborInPlace(Arrays.copyOf(x, x.length), y);
	}
	
	public static NearestNeighbor newNearestNeighborInPlace(double[] x, double[] y) {
		checkXYDimensions(x, y);
		checkMinXLength(x, 2);
		
		NearestNeighbor nh = new NearestNeighbor(x, y);
		nh._calculator = new NeighborCalculator() {
			
			@Override
			public int calculateNeighborIndex(double t) {
				return t >= 0.5 ? 1 : 0;
			}
		};
		return nh;
	}
	
	public static NearestNeighbor newNextNeighbor(double[] x, double[] y) {
		checkXYDimensions(x, y);
		checkMinXLength(x, 2);
		
		NearestNeighbor nh = new NearestNeighbor(Arrays.copyOf(x, x.length), y);
		nh._calculator = new NeighborCalculator() {
			
			@Override
			public int calculateNeighborIndex(double t) {
				return 1;
			}
		};
		return nh;
	}
	
	public static NearestNeighbor newPreviousNeighbor(double[] x, double[] y) {
		checkXYDimensions(x, y);
		checkMinXLength(x, 2);
		
		NearestNeighbor nh = new NearestNeighbor(Arrays.copyOf(x, x.length), y);
		nh._calculator = new NeighborCalculator() {
			
			@Override
			public int calculateNeighborIndex(double t) {
				return 0;
			}
		};
		return nh;
	}

	@Override
	public double evaluateSegmentAt(int index, double x) {
		double t = (x - _x[index]) / (_x[index + 1] - _x[index]);
		return _y[index + _calculator.calculateNeighborIndex(t)];
	}

	@Override
	public UnivariateFunction getSegment(int index) {
		final double yi = this.evaluateSegmentAt(index, _x[index]);
		UnivariateFunction fn = new UnivariateFunction() {
			
			@Override
			public double evaluateAt(double x) {
				return yi;
			}
		};
		return fn;
	}

}
