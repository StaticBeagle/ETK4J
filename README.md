# ETK4J
Engineering Toolkit for Java

This is a project that I work on my spare time. The purpose of the project
is to have a centralized library that's 100% in Java that allows one to prototype solutions
to engineering problems. 

The linear algebra part of the library was taken from [Jama](https://math.nist.gov/javanumerics/jama/)
and the internal data representation was changed from 2d array of doubles to a 1d array of doubles
and the values are accessed via an offset.

Some algorithms in the library are not state-of-the-art, but they should be good enough in terms of accuracy
and speed for many applications.

The library also uses code translated from SciPy and Numpy especially the analog filter part.

There's a set of examples that show how to use some classes contained in the library. The examples
can be found

    src/main/java

# Requirements
---
JDK 1.8+
