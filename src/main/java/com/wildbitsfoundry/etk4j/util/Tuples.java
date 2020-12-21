package com.wildbitsfoundry.etk4j.util;

public final class Tuples {

	private Tuples() {
	}

	public static class Tuple2<T1, T2> {
		private T1 item1;
		private T2 item2;

		public Tuple2(T1 item1, T2 item2) {
			this.item1 = item1;
			this.item2 = item2;
		}

		public T1 getItem1() {
			return item1;
		}

		public T2 getItem2() {
			return item2;
		}

		@Override
		public boolean equals(Object o) {
			if (this == o) return true;
			if (o == null || getClass() != o.getClass()) return false;

			Tuple2<?, ?> tuple2 = (Tuple2<?, ?>) o;

			if (!item1.equals(tuple2.item1)) return false;
			return item2.equals(tuple2.item2);
		}

		@Override
		public int hashCode() {
			int result = item1.hashCode();
			result = 31 * result + item2.hashCode();
			return result;
		}
	}

	public static class Tuple3<T1, T2, T3> extends Tuple2<T1, T2> {
		private T3 item3;

		public Tuple3(T1 item1, T2 item2, T3 item3) {
			super(item1, item2);
			this.item3 = item3;
		}

		public T3 getItem3() {
			return item3;
		}

		@Override
		public boolean equals(Object o) {
			if (this == o) return true;
			if (o == null || getClass() != o.getClass()) return false;
			if (!super.equals(o)) return false;

			Tuple3<?, ?, ?> tuple3 = (Tuple3<?, ?, ?>) o;

			return item3.equals(tuple3.item3);
		}

		@Override
		public int hashCode() {
			int result = super.hashCode();
			result = 31 * result + item3.hashCode();
			return result;
		}
	}
}
