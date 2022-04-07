# ETK4J

## Engineering Toolkit for Java

The purpose of this project is to create a library that can be used to prototype solutions to engineering problems. One 
of the main goals of the library is to be a 100% in Java. Some algorithms in the library are not state-of-the-art, but 
they should be good enough in terms of accuracy and speed for many applications.

The linear algebra part of the library is based on [Jama](https://math.nist.gov/javanumerics/jama/).
The main difference between [Jama](https://math.nist.gov/javanumerics/jama/) and this project is that the internal
representation of the data was changed from a 2d array of doubles to a 1d array of doubles and the values are accessed
using an offset. Other matrix methods were added as well.

The library also uses code translated from [SciPy](https://scipy.org/) and [NumPy](https://numpy.org/). 
Please see [SciPy](https://github.com/StaticBeagle/ETK4J/blob/master/SciPy).

Last but not least, this project includes code that was translated from [numal](https://github.com/JeffBezanson/numal),
and also from [Math.NET](https://www.mathdotnet.com/) please see [Math.NET](https://github.com/StaticBeagle/ETK4J/blob/master/Math.NET.txt).

Some parts of the library do not follow Java naming conventions and the main reason for this is to align with the more familiar Matlab and SciPy syntax.
Arrays are mainly used throughout the library in order to use native doubles but the use of `List`s is encouraged. 

There's a set of examples that show how to use some classes contained in the library. The examples
can be found in:

    src/main/java

# Requirements
---
JDK 1.8+
