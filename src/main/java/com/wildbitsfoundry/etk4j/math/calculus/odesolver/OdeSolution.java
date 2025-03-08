package com.wildbitsfoundry.etk4j.math.calculus.odesolver;

import java.util.ArrayList;
import java.util.List;

public class OdeSolution {

    private final List<DenseOutput> solutionList;
    private final List<Double> t;

    public OdeSolution() {
        solutionList = new ArrayList<>();
        t = new ArrayList<>();
    }

//    public void addDenseOutput(DenseOutput denseOutput) {
//        solutionList.add(denseOutput);
//        t.add(denseOutput.)
//    }
//
//    public int findSegmentIndex(double x) {
//        int index = Arrays.binarySearch(this.x, x);
//        return index < 0.0 ? -(index + 2) : Math.min(index, this.x.length - 2);
//    }
//
//    public final double evaluateAt(double x) {
//        x += 0.0;	// convert -0.0 to 0.0
//        if (x >= x0 && x <= xn) {
//            int index = this.findSegmentIndex(x);
//            return this.evaluateAt(index, x);
//        }
//        return this.extrapolate(x);
//    }

    public final double[] evaluateAt(double[] x) {
        final int n = x.length;
        double[] yi = new double[n];
        for(int i = 0; i < n; ++i) {
            yi[i] = this.evaluateAt(x[i]);
        }
        return yi;
    }

    public double evaluateAt(double t) {
        return 0;
    }
}
