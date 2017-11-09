package com.wildbitsfoundry.etk4j.random;

public class Random {
	private long _seed = System.nanoTime();

	public Random() {
	}

	public long nextLong() {
		long x = _seed;
		x ^= (x << 21);
		x ^= (x >>> 35);
		x ^= (x << 4);
		_seed = x;
		return x;
	}
}
