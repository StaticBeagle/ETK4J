package com.wildbitsfoundry.etk4j.math.interpolation;

public enum ExtrapolationMethod {
	Throw,
	Linear,
	Natural,
	Periodic,
	ClampToNaN,
	ClampToZero,
	ClampToEndPoint
}
