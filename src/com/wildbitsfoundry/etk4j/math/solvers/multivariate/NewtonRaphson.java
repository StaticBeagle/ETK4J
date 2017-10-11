package com.wildbitsfoundry.etk4j.math.solvers.multivariate;

public class NewtonRaphson 
{
	
//	public static double[] Newton(double[] x0, IMultivariateFunction[] functions, IMultivariateFunction[][] jacobian)
//	{
//		int maxiter = 100;
//		double tol = 1e-10;
//		double maxval = 20e4;
//		
//
//		
//		int rows = jacobian.length;
//		int cols = jacobian[0].length;
//		double[] xcurrent = x0;
//		double[] xfinal = new double[cols];
//		double[] functionValues = new double[cols];
//		double[][] jacobianMatrixValues = new double[rows][cols];
//		while(maxiter > 0)
//		{
//			evaluateJacobian(xcurrent, jacobianMatrixValues, jacobian);
//			LU matrix = new LU(jacobianMatrixValues);
//			if(Double.compare(Math.abs(matrix.det()), tol) < 0)
//			{
//				// bail out the matrix is singular
//				return null; //<-- we'll think about this later
//			}
//			
//			// evaluate functions at current x
//			eval(xcurrent, functionValues, functions);
//			double[] residuals = matrix.solve(functionValues);
//			
//			for(int i = 0; i < cols; i++)
//			{
//				xfinal[i] = xcurrent[i] - residuals[i];
//			}
//			eval(xfinal, functionValues, functions);
//			
//			boolean done = true;
//			for(int i = 0; i < cols; i++)
//			{
//				done &= Double.compare(Math.abs(residuals[i]), tol) < 0;
//			}
//			if(done)
//			{
//				
//				return xfinal;
//			}
//			
//			boolean diverged = false;
//			for(int i = 0; i < cols; i++)
//			{
//				diverged |= Double.compare(Math.abs(xcurrent[i]), maxval) > 0;
//			}
//			if(diverged)
//			{
//				return null;
//			}
//			
//			maxiter--;
//			xcurrent = xfinal;
//			
//		}
//		return null;
//	}
//	
//	public static class SolverResults
//	{
//		public boolean Converged;
//		public double[] Results;
//		public int NoIterations;
//		public double Error;
//		public String Status;
//	}
//	
//	public static double[] Newton(final double[] x0, IMultivariateFunction[] functions)
//	{
// 		int maxiter = 100;
//		double tol = 1e-15;
//		double maxval = 20e4;
//		double step = 0.001;
//		double h = 1.0 / step;
//
//		
//		int rows = functions.length;
//		int cols = x0.length;
//		double[] xold = new double[cols];
//		double[] xfinal = new double[cols];
//		double[] residuals = new double[cols];
//		double[] functionValues = new double[cols];
//		double[] functionPrimeValues = new double[cols];
//		double[][] jacobianMatrixValues = new double[rows][cols];
//		
//		System.arraycopy(x0, 0, xold, 0, rows);
//		while(maxiter > 0)
//		{
//			// Compute Jacobian
//			for(int i = 0; i < cols; i++)
//			{
//				xold[i] = xold[i] + step;
//				eval(xold, functionPrimeValues, functions);
//				xold[i] = xold[i] - step;
//				eval(xold, functionValues, functions);
//				for(int j = 0; j < rows; j++)
//				{
//					jacobianMatrixValues[j][i] = (functionPrimeValues[j] - functionValues[j]) * h;
//				}
//			}
//
//			LU matrix = new LU(jacobianMatrixValues);
//			if(Double.compare(Math.abs(matrix.det()), tol) < 0)
//			{
//				// bail out the matrix is singular
//				return null; //<-- we'll think about this later
//			}
//
//			residuals = matrix.solve(functionValues);
//			for(int i = 0; i < cols; i++)
//			{
//				xfinal[i] = xold[i] - residuals[i];
//			}
//			eval(xfinal, functionValues, functions);
//			
//			double relErrorNorm = 0.0;
//			double errorNorm = 0.0;
//			for(int i = 0; i < cols; i++)
//			{
//				// calculate the norm2 of the relative error
//				relErrorNorm += (xfinal[i] - xold[i]) * (xfinal[i] - xold[i]);
//				errorNorm += xfinal[i] * xfinal[i];
//			}
//			double error = Math.sqrt(relErrorNorm) / Math.sqrt(errorNorm);
//			if(Double.compare(error, tol) <= 0)
//			{
//				System.out.println(100 - maxiter);
//				System.out.println(error);
//				return xfinal;
//			}
//			
//			boolean diverged = false;
//			for(int i = 0; i < cols; i++)
//			{
//				diverged |= Double.compare(Math.abs(xold[i]), maxval) > 0;
//			}
//			if(diverged)
//			{
//				return null;
//			}
//			
//			maxiter--;
//			System.arraycopy(xfinal, 0, xold, 0, rows);			
//		}
//		return null;
//	}
//	
//	public static void evaluateJacobian(double[] x, double[][] out, IMultivariateFunction[][] jacobian)
//	{
//		int rows = out.length;
//		int cols = out[0].length;
//		
//		for(int i = 0; i < rows; i++)
//		{
//			for(int j = 0; j < cols; j++)
//			{
//				out[i][j] = jacobian[i][j].EvaluateAt(x);
//			}
//		}
//	}
//	
//	public static void computeJacobian(double[] x, double[][] out, IMultivariateFunction[] functions, double step)
//	{
//		int rows = out.length;
//		int cols = out[0].length;
//		double h = 1.0 / step;
//		double[] fxx = new double[rows];
//		double[] fx = new double[rows];
//		for(int i = 0; i < cols; i++)
//		{
//			x[i] = x[i] + step;
//			eval(x, fxx, functions);
//			x[i] = x[i] - step;
//			eval(x, fx, functions);
//			for(int j = 0; j < rows; j++)
//			{
//				out[j][i] = (fxx[j] - fx[j]) * h;
//			}
//		}
//	}
//	
//	public static void eval(double[]x, double[] out, IMultivariateFunction[] functions)
//	{
//		int length = x.length;
//		for(int i = 0; i < length; i++)
//		{
//			out[i] = functions[i].EvaluateAt(x);
//		}
//	}
//
//	public static void main(String[] args) 
//	{
////		final double vin = 30;
////		final double vb = 0;
////		final double r1 = 10;
////		final double g1 = 1 / r1;
////		final double r2 = 10;
////		final double g2 = 1 / r2;
////		final double L = 1e-3;
////		final double Y = 1 / L;
////		final double tsw = 0.5e-3;
////		
////		// x[0] = ip
////		// x[1] = D
////		// x[2] = iv
////		IMultivariateFunction f1 = new IMultivariateFunction()
////		{
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return -vb * g2 + (x[0] * r2 + vb) * g2 * Math.exp(-r2 * Y * (1 - x[1]) * tsw);
////			}
////		};
////		
////		IMultivariateFunction f2 = new IMultivariateFunction()
////		{
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return (vin - vb) * g1 - (vin - vb - x[2] * r1) * g1 * Math.exp(-r1 * Y * x[1] * tsw); 
////			}
////		};
////		
////		IMultivariateFunction f3 = new IMultivariateFunction()
////		{
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return 1 / tsw * ((vin - vb) * g1 * x[1] * tsw + L * g1 * g1 * (vin - vb - x[2] * r1) * (-1 + Math.exp(-r1 * Y * x[1] * tsw))
////						- g2 * vb * (1 - x[1]) * tsw - L * g2 * g2 *(x[0] * r2 - vb) * (-1 + Math.exp(-r2 * Y * (1 - x[1]) * tsw))); 
////			}
////		};
////		
////		IMultivariateFunction df1x1 = new IMultivariateFunction()
////		{
////
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return -1;
////			}
////		};
////		
////		IMultivariateFunction df1x2 = new IMultivariateFunction()
////		{
////
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return Math.exp(-r2 * Y * (1 - x[1]) * tsw);
////			}
////		};	
////		
////		IMultivariateFunction df1x3 = new IMultivariateFunction()
////		{
////
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return Y * (x[0] * r2 + vb) * tsw * Math.exp(-r1 * Y * x[1] * tsw);
////			}
////		};	
////		
////		IMultivariateFunction df2x1 = new IMultivariateFunction()
////		{
////
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return g1 * Math.exp(-r1 * Y * x[1] * tsw);
////			}
////		};
////		
////		IMultivariateFunction df2x2 = new IMultivariateFunction()
////		{
////
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return -1;
////			}
////		};	
////		
////		IMultivariateFunction df2x3 = new IMultivariateFunction()
////		{
////
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return Y * (vin - x[2] * r1 - vb) * tsw * Math.exp(-r1 * Y * x[1] * tsw);
////			}
////		};	
////		
////		IMultivariateFunction df3x1 = new IMultivariateFunction()
////		{
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return -L * g1 / tsw * (-1 + Math.exp(-r1 * Y * x[1] * tsw));
////			}
////		};
////		
////		IMultivariateFunction df3x2 = new IMultivariateFunction()
////		{
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return -L * g2 / tsw * (-1 + Math.exp(-r2 * Y * (1 - x[1]) * tsw));
////			}
////		};	
////		
////		IMultivariateFunction df3x3 = new IMultivariateFunction()
////		{
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return (vin - vb) * g1 + (vb + x[2] * r1 - vin) * g1 *  Math.exp(-r1 * Y * x[1] * tsw) + (x[0] + vb) * g2 * Math.exp(-r2 * Y * (1 - x[1]) * tsw);
////			}
////		};	
////		IMultivariateFunction[] functions = new IMultivariateFunction[] {f1, f2, f3};
////	
////		IMultivariateFunction[][] jacobian = new IMultivariateFunction[][] {{df1x1, df1x2, df1x3}, {df2x1, df2x2, df2x3}, {df3x1, df3x2, df3x3}};
////		IMultivariateFunction f1 = new IMultivariateFunction()
////		{
////
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return x[0] * x[0] + x[1] * x[1] - 50;
////			}
////		};
////		
////		IMultivariateFunction f2 = new IMultivariateFunction()
////		{
////
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return x[0] * x[1] -25;
////			}
////		};
////		
////		IMultivariateFunction df1x1 = new IMultivariateFunction()
////		{
////
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return 2 * x[0];
////			}
////		};
////		
////		IMultivariateFunction df1x2 = new IMultivariateFunction()
////		{
////
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return 2 * x[1];
////			}
////		};	
////		
////		IMultivariateFunction df2x1 = new IMultivariateFunction()
////		{
////
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return x[1];
////			}
////		};
////		
////		IMultivariateFunction df2x2 = new IMultivariateFunction()
////		{
////
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return x[0];
////			}
////		};	
////		
////		IMultivariateFunction[] functions = new IMultivariateFunction[] {f1, f2};
//	
//		//IMultivariateFunction[][] jacobian = new IMultivariateFunction[][] {{df1x1, df1x2}, {df2x1, df2x2}};
//		IMultivariateFunction f1 = new IMultivariateFunction()
//		{
//
//			@Override
//			public double EvaluateAt(double[] x) 
//			{
//				return x[0] * x[0] - 2 * x[0] +  x[1] * x[1] - x[2] + 1;
//			}
//		};
//		
//		IMultivariateFunction f2 = new IMultivariateFunction()
//		{
//
//			@Override
//			public double EvaluateAt(double[] x) 
//			{
//				return x[0] * x[1] * x[1] - x[0] - 3 * x[1] + x[1] * x[2] + 2;
//			}
//		};
//		
//		IMultivariateFunction f3 = new IMultivariateFunction()
//		{
//
//			@Override
//			public double EvaluateAt(double[] x) 
//			{
//				return x[0] * x[2] * x[2] - 3 * x[2] + x[1] * x[2] * x[2] + x[0] * x[1]; 
//			}
//		};
//
////		IMultivariateFunction f1 = new IMultivariateFunction()
////		{
////
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return 3 * x[0] - Math.cos(x[1] * x[2]) - 3.0 / 2.0;
////			}
////		};
////		
////		IMultivariateFunction f2 = new IMultivariateFunction()
////		{
////
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return 4 * x[0] * x[0] - 625 * x[1] * x[1] + 2 * x[2] - 1;
////			}
////		};
////		
////		IMultivariateFunction f3 = new IMultivariateFunction()
////		{
////
////			@Override
////			public double EvaluateAt(double[] x) 
////			{
////				return 20 * x[2] + Math.exp(-x[0] * x[1]) + 9;
////			}
////		};
//		
//		
//		IMultivariateFunction[] functions = new IMultivariateFunction[] {f1, f2, f3};
//		double[] initialguess = new double[] {1d,2d,3d};
//		double[] solution = Newton(initialguess, functions);
//
//		for(double val : solution)
//		{
//			System.out.println(val);
//		}
//		System.out.println("done");
//
//	}

}
