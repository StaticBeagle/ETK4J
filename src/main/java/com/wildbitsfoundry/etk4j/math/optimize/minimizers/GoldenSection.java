package com.wildbitsfoundry.etk4j.math.optimize.minimizers;

import java.util.function.BiFunction;

public final class GoldenSection {

    private GoldenSection() {}

    // TODO
    // Maybe improve it and roll out our own implementation
    // https://www.mathworks.com/matlabcentral/fileexchange/25919-golden-section-method-algorithm
    public static double goldenSectionMinimizer(BiFunction<Double, Object[], Double> func,
                                                double a, double b, double tol, int maxIter, Object... params) {
        double gold = (Math.sqrt(5.0) - 1.0) / 2.0;

        double x1 = a + (1 - gold) * (b - a);
        double x2 = a + gold * (b - a);

        double fx1 = func.apply(x1, params);
        double fx2 = func.apply(x2, params);

        int k = 1;
        while ((Math.abs(b - a) > tol) && (k < maxIter)) {
            k = k + 1;
            if (fx1 < fx2) {
                b = x2;
                x2 = x1;
                x1 = a + (1 - gold) * (b - a);
            } else {
                a = x1;
                x1 = x2;
                x2 = a + gold * (b - a);
            }
            fx1 = func.apply(x1, params);
            fx2 = func.apply(x2, params);
        }
        if (fx1 < fx2)
            return x1;
        return x2;
    }
}
