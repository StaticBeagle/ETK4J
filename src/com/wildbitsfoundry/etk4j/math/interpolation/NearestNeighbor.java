package com.wildbitsfoundry.etk4j.math.interpolation;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.functions.PiecewiseFunction;
import com.wildbitsfoundry.etk4j.util.ArrayUtils;

public class NearestNeighbor extends PiecewiseFunction {
	
	interface NeighborCalculator {
		int calculateNeighborIndex(double t);
	}
	
	private double[] _y = null;
	private NeighborCalculator _calculator = null;
	protected NearestNeighbor(double[] x, double[] y) {
		super(x, y[0], y[y.length - 1]);
		_y = Arrays.copyOf(y, y.length);
		if(!ArrayUtils.isAscending(x)) {
			throw new IllegalArgumentException("x must be monotonically increasing");
		}
		this.setExtrapolationMethod(ExtrapolationMethod.ClampToEndPoint);
	}
	
	public static NearestNeighbor newNearestNeighbor(double[] x, double[] y) {
		return newNearestNeighborInPlace(Arrays.copyOf(x, x.length), y);
	}
	
	public static NearestNeighbor newNearestNeighborInPlace(double[] x, double[] y) {
		Splines.checkXYDimensions(x, y);
		Splines.checkMinkXLength(x, 2);
		
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
		Splines.checkXYDimensions(x, y);
		Splines.checkMinkXLength(x, 2);
		
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
		Splines.checkXYDimensions(x, y);
		Splines.checkMinkXLength(x, 2);
		
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
	protected double getValueAt(int i, double x) {
		double t = (x - _x[i]) / (_x[i + 1] - _x[i]);
		return _y[i + _calculator.calculateNeighborIndex(t)];
	}

}
