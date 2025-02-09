package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertTrue;

public class ComplexCholeskyDenseTest {

    @Test
    public void testGetR() {
        ComplexMatrixDense Arg = ComplexMatrixDense.from2DArray(new Complex[][]{
                {new Complex(0.9767407521087846, 0.4711470863081668), new Complex(0.2794308404798361, 0.4587985898568102), new Complex(0.1685348957730604, 0.2220415588608931), new Complex(0.8646965674039528, 0.5496275298999502)},
                {new Complex(0.7198733467511716, 0.9814006625666408), new Complex(0.01506799866339292, 0.666774973911552), new Complex(0.750545968546426, 0.4842583141948102), new Complex(0.8173900466200722, 0.9463910637201288)},
                {new Complex(0.3444764277139777, 0.7365028309538496), new Complex(0.8429046576872399, 0.1345623113997281), new Complex(0.6191817846917982, 0.174500592061623), new Complex(0.1931980127831494, 0.3833402493069025)},
                {new Complex(0.2529620785541479, 0.8308769763271352), new Complex(0.7385598627121508, 0.850513285400206), new Complex(0.3669833799842998, 0.3169573549259639), new Complex(0.4630012679866645, 0.4444797942481419)}
        });

        // Making it symmetric
        Arg.multiplyEquals(Arg.conjugateTranspose());
        ComplexCholeskyDecompositionDense choleskyDecompositionDense = new ComplexCholeskyDecompositionDense(Arg);

        Complex[] expected = {new Complex(1.6099928119998814, 0), new Complex(1.8239927308960766, -0.6725838511782638),
                new Complex(0.9326675288675592, -0.20196491214994378), new Complex(1.2497046390097175, -0.4304200565982681),
                new Complex(0, 0), new Complex(0.7130709548221029, -0),
                new Complex(0.4293064714522471, 0.20662989447609814), new Complex(0.3293349304457535, 0.12482789775498358),
                new Complex(0, 0), new Complex(0, 0),
                new Complex(0.9220352270895488, 0), new Complex(0.6060228453857976, -0.6447987779995583),
                new Complex(0, 0), new Complex(0, 0),
                new Complex(0, 0), new Complex(0.12714301598503266, 0)};
        assertArrayEquals(expected, choleskyDecompositionDense.getR().getArray());
    }

    @Test
    public void testGetOriginalMatrixBack() {
        ComplexMatrixDense Arg = ComplexMatrixDense.from2DArray(new Complex[][]{
                {new Complex(0.9767407521087846, 0.4711470863081668), new Complex(0.2794308404798361, 0.4587985898568102), new Complex(0.1685348957730604, 0.2220415588608931), new Complex(0.8646965674039528, 0.5496275298999502)},
                {new Complex(0.7198733467511716, 0.9814006625666408), new Complex(0.01506799866339292, 0.666774973911552), new Complex(0.750545968546426, 0.4842583141948102), new Complex(0.8173900466200722, 0.9463910637201288)},
                {new Complex(0.3444764277139777, 0.7365028309538496), new Complex(0.8429046576872399, 0.1345623113997281), new Complex(0.6191817846917982, 0.174500592061623), new Complex(0.1931980127831494, 0.3833402493069025)},
                {new Complex(0.2529620785541479, 0.8308769763271352), new Complex(0.7385598627121508, 0.850513285400206), new Complex(0.3669833799842998, 0.3169573549259639), new Complex(0.4630012679866645, 0.4444797942481419)}
        });

        // Making it symmetric
        Arg.multiplyEquals(Arg.conjugateTranspose());
        ComplexCholeskyDecompositionDense choleskyDecompositionDense = new ComplexCholeskyDecompositionDense(Arg);

        ComplexMatrixDense R = choleskyDecompositionDense.getR();
        ComplexMatrixDense A = R.conjugateTranspose().multiply(R);
        double[] argMag = Arrays.stream(Arg.getArray()).mapToDouble(Complex::abs).toArray();
        double[] aMag = Arrays.stream(A.getArray()).mapToDouble(Complex::abs).toArray();
        assertArrayEquals(argMag, aMag, 1e-12);
    }

    @Test
    public void testIsSPD() {
        ComplexMatrixDense Arg = ComplexMatrixDense.from2DArray(new Complex[][]{
                {new Complex(0.9767407521087846, 0.4711470863081668), new Complex(0.2794308404798361, 0.4587985898568102), new Complex(0.1685348957730604, 0.2220415588608931), new Complex(0.8646965674039528, 0.5496275298999502)},
                {new Complex(0.7198733467511716, 0.9814006625666408), new Complex(0.01506799866339292, 0.666774973911552), new Complex(0.750545968546426, 0.4842583141948102), new Complex(0.8173900466200722, 0.9463910637201288)},
                {new Complex(0.3444764277139777, 0.7365028309538496), new Complex(0.8429046576872399, 0.1345623113997281), new Complex(0.6191817846917982, 0.174500592061623), new Complex(0.1931980127831494, 0.3833402493069025)},
                {new Complex(0.2529620785541479, 0.8308769763271352), new Complex(0.7385598627121508, 0.850513285400206), new Complex(0.3669833799842998, 0.3169573549259639), new Complex(0.4630012679866645, 0.4444797942481419)}
        });

        // Making it symmetric
        Arg.multiplyEquals(Arg.conjugateTranspose());
        ComplexCholeskyDecompositionDense choleskyDecompositionDense = new ComplexCholeskyDecompositionDense(Arg);

       assertTrue(choleskyDecompositionDense.isSPD());
    }


    @Test
    public void testSolve() {
        ComplexMatrixDense Arg = ComplexMatrixDense.from2DArray(new Complex[][]{
                {new Complex(0.9767407521087846, 0.4711470863081668), new Complex(0.2794308404798361, 0.4587985898568102), new Complex(0.1685348957730604, 0.2220415588608931), new Complex(0.8646965674039528, 0.5496275298999502)},
                {new Complex(0.7198733467511716, 0.9814006625666408), new Complex(0.01506799866339292, 0.666774973911552), new Complex(0.750545968546426, 0.4842583141948102), new Complex(0.8173900466200722, 0.9463910637201288)},
                {new Complex(0.3444764277139777, 0.7365028309538496), new Complex(0.8429046576872399, 0.1345623113997281), new Complex(0.6191817846917982, 0.174500592061623), new Complex(0.1931980127831494, 0.3833402493069025)},
                {new Complex(0.2529620785541479, 0.8308769763271352), new Complex(0.7385598627121508, 0.850513285400206), new Complex(0.3669833799842998, 0.3169573549259639), new Complex(0.4630012679866645, 0.4444797942481419)}
        });

        // Making it symmetric
        Arg.multiplyEquals(Arg.conjugateTranspose());

        ComplexCholeskyDecompositionDense choleskyDecompositionDense = new ComplexCholeskyDecompositionDense(Arg);
        ComplexMatrixDense b = ComplexMatrixDense.from2DArray(new Complex[][] {
                {new Complex(90, 90)},
                {new Complex(80, 80)},
                {new Complex(70, 70)},
                {new Complex(60, 60)}
        });

        ComplexMatrixDense X = choleskyDecompositionDense.solve(b);
        Complex[] expected = {new Complex(403.61865072620583, 2427.3094092712613), new Complex(-1276.4687915825514, -1323.8099753278434),
                new Complex(1602.853981799121, 3633.4514965471276), new Complex(1604.0632231016468, -3737.3817702615556)};
        assertArrayEquals(expected, X.getArray());
    }

    @Test
    public void testBalance() {
        ComplexMatrixDense Arg = ComplexMatrixDense.from2DArray(new Complex[][]{
                {new Complex(1, 50), new Complex(100, 150), new Complex(10000, 50000)},
                {new Complex(0.1, 0.1), new Complex(1, 1), new Complex(100, 500)},
                {new Complex(0.0001, 0.0001), new Complex(0.1, 0.5), new Complex(1, 15)}
        });

        Complex[] expected = new Complex[] {new Complex(1, 50), new Complex(1.5625, 2.34375), new Complex(2.44140625, 12.20703125),
        new Complex(6.4, 6.4), new Complex(1, 1), new Complex(1.5625, 7.8125),
        new Complex(0.4096, 0.4096), new Complex(6.4, 32), new Complex(1, 15)};

        assertArrayEquals(expected, Arg.balance().getArray());
    }
}
