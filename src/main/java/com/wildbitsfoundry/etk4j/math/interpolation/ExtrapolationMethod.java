package com.wildbitsfoundry.etk4j.math.interpolation;

// TODO this all should be all capitalized
public enum ExtrapolationMethod {
	THROW,
	LINEAR,
	NATURAL,
	PERIODIC,
	CLAMP_TO_NAN,
	CLAMP_TO_ZERO,
	CLAMP_TO_END_POINT
}
