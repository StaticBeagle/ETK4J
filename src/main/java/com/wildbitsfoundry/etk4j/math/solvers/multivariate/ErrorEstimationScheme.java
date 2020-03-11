package com.wildbitsfoundry.etk4j.math.solvers.multivariate;

import com.wildbitsfoundry.etk4j.util.NumArrays;

public enum ErrorEstimationScheme {
	// Norm 1
	SUM_ABS_ERROR {				
		@Override
		double calculateRelativeError(double[] x, double[] y) {
			return NumArrays.norm1(NumArrays.subtract(x, y));
		}

		@Override
		double calculateMaxNormOfError(double[] x, double[] y) {
			return Math.min(NumArrays.norm1(x), NumArrays.norm1(y));
		}
	},
	// Infinite Norm
	MAX_ABS_ERROR {				
		@Override
		double calculateRelativeError(double[] x, double[] y) {
			return NumArrays.normInf(NumArrays.subtract(x, y));
		}

		@Override
		double calculateMaxNormOfError(double[] x, double[] y) {
			return Math.min(NumArrays.normInf(x), NumArrays.normInf(y));
		}
	},
	// Norm 2
	SQRT_SUM_ABS_ERROR_SQUARED {	
		@Override
		double calculateRelativeError(double[] x, double[] y) {
			return NumArrays.distance(x, y);
		}

		@Override
		double calculateMaxNormOfError(double[] x, double[] y) {
			return Math.min(NumArrays.norm2(x), NumArrays.norm2(y));
		}
	};

	abstract double calculateRelativeError(double[] x, double[] y);
	abstract double calculateMaxNormOfError(double[] x, double[] y);
}