package com.wildbitsfoundry.etk4j.math.solvers.multivariate;

public enum ErrorEstimationScheme {
	SUM_ABS_ERROR,				// Norm 1
	MAX_ABS_ERROR,				// Infinite Norm
	SQRT_SUM_ABS_ERROR_SQUARED,	// Norm 2
}
