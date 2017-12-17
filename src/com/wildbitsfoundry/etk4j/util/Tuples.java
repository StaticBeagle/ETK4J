package com.wildbitsfoundry.etk4j.util;

public final class Tuples {
	private Tuples() {
	}

	public static class Tuple2<T1, T2> {
		public T1 Item1;
		public T2 Item2;
		
		public Tuple2(T1 item1, T2 item2) {
			Item1 = item1;
			Item2 = item2;
		}
	}
	
	public static class Tuple3<T1, T2, T3> extends Tuple2<T1, T2> {
		public T3 Item3;
		
		public Tuple3(T1 item1, T2 item2, T3 item3) {
			super(item1, item2);
			Item3 = item3;
		}
	}
}
