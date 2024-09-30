package com.wildbitsfoundry.etk4j.util.data;

import org.junit.Test;

import java.util.Arrays;
import java.util.stream.IntStream;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

public class DynamicArrayTest {

    @Test
    public void testDefaultSettings() {
        DynamicArray dynamicArray = new DynamicArray();
        assertEquals(0, dynamicArray.size());
    }

    @Test
    public void testArrayGrowsInSize() {
        DynamicArray dynamicArray = new DynamicArray();
        double[] data = IntStream.range(0, 20).asDoubleStream().toArray();
        Arrays.stream(data).forEach(dynamicArray::add);
        assertArrayEquals(data, IntStream.range(0, 20).mapToDouble(dynamicArray::get).toArray(), 1e-12);
    }
}
