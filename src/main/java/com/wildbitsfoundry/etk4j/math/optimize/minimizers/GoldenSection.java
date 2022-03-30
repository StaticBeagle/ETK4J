package com.wildbitsfoundry.etk4j.math.optimize.minimizers;

import java.util.function.BiFunction;

public final class GoldenSection {

    private final BiFunction<Double, Object[], Double> function;
    private final double a;
    private final double b;
    private final Object[] params;
    private double tol = 1e-5;
    private int maxNumberOfIterations = 500;

    public GoldenSection(BiFunction<Double, Object[], Double> function, double a, double b, Object... params) {
        this.function = function;
        this.a = a;
        this.b = b;
        this.params = params;
    }

    public GoldenSection tolerance(double tol) {
        this.tol = tol;
        return this;
    }

    public GoldenSection iterationLimit(int limit) {
        this.maxNumberOfIterations = limit;
        return this;
    }

    // https://www.mathworks.com/matlabcentral/fileexchange/25919-golden-section-method-algorithm
    public MinimizerResults<Double> minimize() {
        double gold = (Math.sqrt(5.0) - 1.0) / 2.0;
        double x0 = this.a;
        double xn = this.b;
        double x1 = x0 + (1 - gold) * (xn - x0);
        double x2 = x0 + gold * (xn - x0);

        double fx1 = function.apply(x1, params);
        double fx2 = function.apply(x2, params);

        int k = 1;
        while ((Math.abs(xn - x0) > tol) && (k < maxNumberOfIterations)) {
            k = k + 1;
            if (fx1 < fx2) {
                xn = x2;
                x2 = x1;
                x1 = x0 + (1 - gold) * (xn - x0);
            } else {
                x0 = x1;
                x1 = x2;
                x2 = x0 + gold * (xn - x0);
            }
            fx1 = function.apply(x1, params);
            fx2 = function.apply(x2, params);
        }
        if(k >= maxNumberOfIterations) {
            if(fx1 < fx2) {
                MinimizerResults<Double> minResults = new MinimizerResults<>();
                minResults.setMinimizerStatus("Maximum number of iterations exceeded");
                minResults.setHasConverged(false);
                minResults.setValue(x1);
                minResults.setFunctionValue(fx1);
                minResults.setNumberOfIterations(k);
                return minResults;
            }
            MinimizerResults<Double> minResults = new MinimizerResults<>();
            minResults.setMinimizerStatus("Maximum number of iterations exceeded");
            minResults.setHasConverged(false);
            minResults.setValue(x2);
            minResults.setFunctionValue(fx2);
            minResults.setNumberOfIterations(k);
            return minResults;
        }
        if (fx1 < fx2) {
            MinimizerResults<Double> minResults = new MinimizerResults<>();
            minResults.setMinimizerStatus("Converged");
            minResults.setHasConverged(true);
            minResults.setValue(x1);
            minResults.setFunctionValue(fx1);
            minResults.setNumberOfIterations(k);
            return minResults;
        }
        MinimizerResults<Double> minResults = new MinimizerResults<>();
        minResults.setMinimizerStatus("Converged");
        minResults.setHasConverged(true);
        minResults.setValue(x2);
        minResults.setFunctionValue(fx2);
        minResults.setNumberOfIterations(k);
        return minResults;
    }
}
