package com.wildbitsfoundry.etk4j.math.solvers.multivariate;

import com.wildbitsfoundry.etk4j.util.DoubleArrays;

public enum ErrorEstimationScheme {
	// Norm 1
	SUM_ABS_ERROR {				
		@Override
		double calculateRelativeError(double[] x, double[] y) {
			return DoubleArrays.norm1(DoubleArrays.subtract(x, y));
		}

		@Override
		double calculateMaxNormOfError(double[] x, double[] y) {
			return Math.min(DoubleArrays.norm1(x), DoubleArrays.norm1(y));
		}
	},
	// Infinite Norm
	MAX_ABS_ERROR {				
		@Override
		double calculateRelativeError(double[] x, double[] y) {
			return DoubleArrays.normInf(DoubleArrays.subtract(x, y));
		}

		@Override
		double calculateMaxNormOfError(double[] x, double[] y) {
			return Math.min(DoubleArrays.normInf(x), DoubleArrays.normInf(y));
		}
	},
	// Norm 2
	SQRT_SUM_ABS_ERROR_SQUARED {	
		@Override
		double calculateRelativeError(double[] x, double[] y) {
			return DoubleArrays.distance(x, y);
		}

		@Override
		double calculateMaxNormOfError(double[] x, double[] y) {
			return Math.min(DoubleArrays.norm2(x), DoubleArrays.norm2(y));
		}
	};

	abstract double calculateRelativeError(double[] x, double[] y);
	abstract double calculateMaxNormOfError(double[] x, double[] y);
}