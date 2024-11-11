package com.wildbitsfoundry.etk4j.util.data;

import org.junit.Test;

import java.util.Arrays;
import java.util.stream.IntStream;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

public class DoubleDynamicArrayTest {

    @Test
    public void testDefaultSettings() {
        DoubleDynamicArray doubleDynamicArray = new DoubleDynamicArray();
        assertEquals(0, doubleDynamicArray.size());
    }

    @Test
    public void testArrayGrowsInSize() {
        DoubleDynamicArray doubleDynamicArray = new DoubleDynamicArray();
        double[] data = IntStream.range(0, 20).asDoubleStream().toArray();
        Arrays.stream(data).forEach(doubleDynamicArray::add);
        assertArrayEquals(data, IntStream.range(0, 20).mapToDouble(doubleDynamicArray::get).toArray(), 1e-12);
    }
}
