// 
/*
Author:		Brian Pendleton
UNCID:		800458057
Date:		10/03/2010
Class:		MATH-6204
Semester:	FALL-2010

Purpose:
The following program is designed to apply a central-difference approximation
to a function f(x) = e^x, with a decreasing sequence of h.  It is calculated
two ways, one using floats and one using doubles.

We will use templates so we don't rewrite the same code twice.  A template is a 
function that can execute for any type, without that type being hard coded.  Example:

int Add(int a,int b) { return a+b;}			// function Without C++ template
float Add(float a, float b) { return a+b;}	// function Without C++ template

These functions can be combined using a template which works for both int and float
template <class T>
T Add(T a, T b)		//C++ function template sample
{
	return a+b;
}

*/
//
#if defined(__WIN32__) || defined(_WIN32_) || defined(__WIN32) || defined(_WIN32) || defined(WIN32)
#include "stdafx.h"
#endif
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

using namespace std;

template <class T>
T fPrime(T, T);

template <class T>
T f(T);

template <class T>
T h(T);

template <class T>
void calcFPrimeAndPrint(T,T);

int main(int argc, char* argv[])
{
	float f_h;
	float f_x = 0.0;
	double d_h;
	double d_x = 0.0;

	cout << "Results using FLOAT" << endl;
	cout << left << setw(5) << "i" << left << setw(20) << "h" << left << setw(20) << "f'(0)" << endl;
	cout << "----------------------------------------------" << endl;

	for ( int i=0; i<=25; i++ )
	{
		f_h = h( (float)i);
		cout << left << setw(5) << i;

		calcFPrimeAndPrint(f_x, f_h);
	}
	cout << endl;
	cout << "Results using DOUBLE" << endl;
	cout << left << setw(5) << "i" << left << setw(20) << "h" << left << setw(20) << "f'(0)" << endl;
	cout << "----------------------------------------------" << endl;
	for ( int i=0; i<=25; i++ )
	{
		d_h = h( (double)i);
		cout << left << setw(5) << i;

		calcFPrimeAndPrint(d_x, d_h);
	}


}

/*
Wrapper function to move this logic out of the main function.
This simply calls the fPrime function with the templated values of x and h.
*/
template <class T>
void calcFPrimeAndPrint(T x, T h)
{
	cout << fixed;
	cout.width(20);
	cout << left << setprecision(8) << h << left << setprecision(8) << fPrime(x, h) << endl;

}


/*
Template class to calculate values of h in the problem.
*/
template <class T>
T h(T i)
{
	return (1 / ((T)2*pow((T)2,i)));
}


/*
Templated class that calculates f'(x) given by the formula stated in the problem.
Because it is templated, it can be calculated for any type 
T (int, float, double, long double, etc...)
*/
template <class T>
T fPrime(T x, T h)
{
	return ( ( f(x+h) - f(x-h) ) / ((T)2*h) );
}


/*
Returns e^x for the Template (double, float, long double, etc....).
*/
template <class T>
T f(T x)
{
	return (exp(x));
}

