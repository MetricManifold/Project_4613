#pragma once


/*
 * this namespace contains functions that will return the value of the corresponding
 * inhomogeneous function for the heat equation
 * in order to set up the constants, one call at the beginning of the program needs to be made
 * that will set up the static constants
 */
namespace inh
{
	/*
	 * returns a constant value
	 * parameter 'c' is given in the first call to this function and is the value
	 * that will be returned on every subsequent call
	 */
	inline double constant(double c = 0)
	{
		static double _c = c;
		return _c;
	}

	/* 
	 * returns a quadratic function with coefficients for both x and y directions
	 * parameters 'x' and 'y' are the coordinates of the node at which to compute this function
	 * parameters 'a' and 'b' are called at the first call to this function and set the constants
	 * used to compute the output on every subsequent call
	 */
	inline double quadratic(double x, double y, double a = 1, double b = 1)
	{
		static double _a = a;
		static double _b = b;

		return _a * x * x + _b * y * y;
	}

}