package com.wildbitsfoundry.etk4j.control;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;
import com.wildbitsfoundry.etk4j.util.NumArrays;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static com.wildbitsfoundry.etk4j.util.ComplexArrays.deepCopy;

public class ZeroPoleGain extends LinearTimeInvariantSystem {
    private final Complex[] zeros;
    private final Complex[] poles;
    private final double gain;

    public ZeroPoleGain(Complex[] zeros, Complex[] poles, double gain) {
        this.zeros = deepCopy(zeros);
        this.poles = deepCopy(poles);
        this.gain = gain;
    }

    public Complex[] getZeros() {
        return deepCopy(zeros);
    }

    public Complex[] getPoles() {
        return deepCopy(poles);
    }

    public double getGain() {
        return gain;
    }

    @Override
    public StateSpace toStateSpace() {
        return new TransferFunction(this).toStateSpace();
    }

    @Override
    public TransferFunction toTransferFunction() {
        return new TransferFunction(this);
    }

    @Override // TODO this should be public to comply with Liskov
    public ZeroPoleGain toZeroPoleGain() {
        return this;
    }

    private static class ComplexRealResults {
        private double[] real;
        private Complex[] complex;

        ComplexRealResults(Complex[] complex, double[] real) {
            this.complex = complex;
            this.real = real;
        }

        public Complex[] getComplex() {
            return complex;
        }

        public double[] getReal() {
            return real;
        }
    }

    private static ComplexRealResults complexReal(Complex[] z) {
        return complexReal(z, 100 * ConstantsETK.DOUBLE_EPS);
    }

    private static ComplexRealResults complexReal(Complex[] z, double tol) {
        if (z.length == 0) {
            return new ComplexRealResults(new Complex[0], new double[0]);
        }

        Arrays.sort(z, Comparator.comparingDouble(Complex::real).thenComparingDouble(o -> Math.abs(o.imag())));

        // Split reals from conjugate pairs
        List<Boolean> zrl = new ArrayList<>();
        for (int i = 0; i < z.length; ++i) {
            if (Math.abs(z[i].imag()) <= tol * z[i].abs()) {
                zrl.add(true);
            } else {
                zrl.add(false);
            }
        }
        double[] zr = new double[(int) zrl.stream().filter(Boolean::booleanValue).count()];
        for (int i = 0, j = 0; i < z.length; ++i) {
            if (zrl.get(i)) {
                zr[j] = z[i].real();
                ++j;
            }
        }

        // Input is entirely real
        if (zr.length == z.length) {
            return new ComplexRealResults(new Complex[0], zr);
        }

        Complex[] zz = new Complex[z.length - zr.length];
        for (int i = 0, j = 0; i < z.length; ++i) {
            if (!zrl.get(i)) {
                zz[j] = z[i];
                ++j;
            }
        }

        List<Complex> zp = new ArrayList<>(zz.length);
        List<Complex> zn = new ArrayList<>(zz.length);
        for (Complex c : zz) {
            if (c.imag() > 0) {
                zp.add(c);
            } else if (c.imag() < 0) {
                zn.add(c);
            }
        }

        if (zp.size() != zn.size()) {
            // TODO throw exception complex values with no matching conjugate
        }

        double[] diffReal = NumArrays.difference(ComplexArrays.real(zp.toArray(new Complex[zp.size()])));
        double[] sameReal = new double[diffReal.length];
        for (int i = 0; i < diffReal.length; ++i) {
            if (diffReal[i] <= tol * zp.get(i).abs()) {
                sameReal[i] = 1.0;
            }
        }

        double[] diffs = NumArrays.difference(NumArrays.concatenateAll(new double[]{0.0}, sameReal, new double[]{0.0}));

        int runStarts = -1;
        for (int i = 0; i < diffs.length; ++i) {
            if (diffs[i] > 0) {
                runStarts = i;
                break;
            }
        }
        int runStops = -1;
        for (int i = 0; i < diffs.length; ++i) {
            if (diffs[i] < 0) {
                runStops = i;
                break;
            }
        }

        if (runStarts >= 0) {
            int start = runStarts;
            int stop = runStops + 1;

            zp.subList(start, stop).sort(Comparator.comparingDouble(o -> Math.abs(o.imag())));
            zn.subList(start, stop).sort(Comparator.comparingDouble(o -> Math.abs(o.imag())));
        }

//            # Check that negatives match positives
//		if any(abs(zp - zn.conj()) > tol * abs(zn)):
//		raise ValueError('Array contains complex value with no matching '
//		'conjugate.')
        for (int i = 0; i < zp.size(); ++i) {
            if (zp.get(i).subtract(zn.get(i).conj()).abs() > tol * zn.get(i).abs()) {
                // throw no matching conjugate
                System.out.println("here");
            }
        }

        Complex[] zc = new Complex[zp.size()];
        // Average out numerical inaccuracy in real vs imag parts of pairs
        for (int i = 0; i < zp.size(); ++i) {
            zc[i] = (zp.get(i).add(zn.get(i).conj())).multiply(0.5);
        }
        return new ComplexRealResults(zc, zr);
    }

    // TODO create enum for realOrComplex
    private int nearestRealComplexIndex(Complex[] fro, Complex to, String which) {
////        """Get the next closest real or complex element based on distance"""
//        order = np.argsort(np.abs(fro - to))
//        mask = np.isreal(fro[order])
//        if which == 'complex':
//        mask = ~mask
//        return order[np.nonzero(mask)[0][0]]
        int[] order = NumArrays.argSort(ComplexArrays.abs(ComplexArrays.subtract(fro, to)));
        // Sort mask array based on order array
        int[] mask = new int[order.length];
        for (int i = 0; i < order.length; ++i) {
            mask[i] = fro[order[i]].isReal() ? 1 : 0;
        }
        if (which.equals("complex")) {
            // mask = ~mask
            mask = Arrays.stream(mask).map(i -> i == 1 ? 0 : 1).toArray();
        }

        // get nonzero indexes
        List<Integer> nonZero = new ArrayList<>();
        for (int i = 0; i < mask.length; ++i) {
            if (mask[i] != 0) {
                nonZero.add(i);
            }
        }
        return order[nonZero.get(0)];

    }

    public static void main(String[] args) {
        Complex[] a = new Complex[]{
                Complex.fromReal(4.0),
                Complex.fromReal(3.0),
                Complex.fromReal(1.0),
                new Complex(2, -2),
                new Complex(2, 2),
                new Complex(2, -1),
                new Complex(2, 1),
                new Complex(2, -1),
                new Complex(2, 1),
                new Complex(1, 1),
                new Complex(1, -1)
        };

        //ComplexRealResults cplx = complexReal(a);

//        Complex[] z1 = new Complex[]{
//                Complex.fromReal(-1.0),
//                new Complex(-0.5, -0.5),
//                new Complex(-0.5, 0.5),
//                Complex.fromReal(4.0),
//                Complex.fromReal(3.0),
//                Complex.fromReal(1.0),
//                new Complex(2, -2),
//                new Complex(2, 2),
//                new Complex(2, -1),
//                new Complex(2, 1),
//                new Complex(1, 1),
//                new Complex(1, -1)
//        };
//        Complex[] p1 = new Complex[]{
//                Complex.fromReal(0.75),
//                new Complex(0.8, 0.1),
//                new Complex(0.8, -0.1),
//                Complex.fromReal(8.0),
//                Complex.fromReal(6.0),
//                Complex.fromReal(2.0),
//                new Complex(4, -4),
//                new Complex(4, 4),
//                new Complex(4, -2),
//                new Complex(4, 2),
//                new Complex(2, 2),
//                new Complex(2, -2)
//        };

        Complex[] z1 = new Complex[0];
        Complex[] p1 = {Complex.fromReal(-1)};

        ZeroPoleGain zpk = new ZeroPoleGain(z1, p1, 1.0);
    }
}
