/* 
HW2.cpp
Brian Pendleton
ID: 800458057

IDE: Visual Studio 2010
Compiled & Executed using Cygwin GNU g++

Description:
HW2 -- Calculate option prices for Euro Call and Put.
*/

#if defined(__WIN32__) || defined(_WIN32_) || defined(__WIN32) || defined(_WIN32) || defined(WIN32)
#include "stdafx.h"
#endif
#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;

/* Here I define some constants to use in the estimation of F(X) */
#if !defined(MATH_CONSTANTS)
#define MATH_CONSTANTS
#define PI 3.14159265358979323846
#define A1 0.319381530
#define A2 -0.356563782
#define A3 1.781477937
#define A4 -1.821255978
#define A5 1.330274429
#endif

double F(double x);
double f(double x);
double EuropeanCall(double S, double K, double r, double delta, double sigma, double T, double t);
double EuropeanPut(double S, double K, double r, double delta, double sigma, double T, double t);
double d1(double S, double K, double r, double delta, double sigma, double T, double t);
double d2(double d1, double sigma, double T, double t);

int main(int argc, char* argv[])
{
	/* Here we declare the initial paramters */
	double K = 100;
	double S_0 = 100;
	double t = 0;
	double T = 1;
	double delta = 0.03;
	double r = 0.025;
	double sigma = 0.0;


	/*  Print the headers for some formatted output */
	cout << "___________ Initial Value of Parameters ___________" << endl << endl;
	cout << setw(10) << "K = " << K << endl;
	cout << setw(10) << "S_0 = " << S_0 << endl;
	cout << setw(10) << "T = " << T << endl;
	cout << setw(10) << "t = " << t << endl;
	cout << setw(10) << "delta = " << delta << endl;
	cout << setw(10) << "r = " << r << "\n\n\n";
	

	/* This prints the column headers for the 10 simulations */
	cout << "_______ Simulated Results of Option Prices _______" << endl << endl;
	cout << setw(10) << "Sigma" << setw(20) << "European Call" << setw(20) << "European Put" << endl;

	/* Start at sigma = .1, and run the simulation. Loop until you've completed sigma = 1.00 */
	for ( sigma = .1; sigma <= 1.00; sigma += .1 )
	{
		/* Call the functions to calculate the price of the Call and Put */
		double V_C = EuropeanCall(S_0, K, r, delta, sigma, T, t);
		double V_P = EuropeanPut(S_0, K, r, delta, sigma, T, t);
		
		cout << fixed;
		cout << setw(10) << setprecision(2) << sigma;
		cout << setw(20) << setprecision(8) << V_C;
		cout << setw(20) << setprecision(8) << V_P << endl;
	}
}

/* This function calculates the value f(x) from the density function */ 
double f(double x)
{
	return ((exp(-(pow(x,2))/2))/sqrt(2*PI));
}

/* This function estimates the Normal Cumulative Distribution.  Note
that if the value of X is negative, we recursively call the function
again using the property F(-x) = 1 - F(x) */
double F(double x)
{
	if ( x < 0 )
	{
		return ( 1 - F(fabs(x)) );
	}
	else 
	{
		double z = ( 1 / (1+.2316419*x));
		return ( 1 - f(x)*z*((((A5*z + A4)*z + A3)*z + A2)*z + A1));  // Uses the constants from #define above.
	}
	
}

/* Calculates the Price of an option.  First invoking the functions to compute d1 and d2.  Then returns
the calculation given in the homework description. */
double EuropeanCall(double S, double K, double r, double delta, double sigma, double T, double t)
{
	double _d1 = d1(S, K, r, delta, sigma, T, t);
	double _d2 = d2(_d1, sigma, T, t);
	return ( S*exp(-delta*(T-t))*F(_d1) - K*exp(-r*(T-t))*F(_d2) );
}

/* Calculates the Price of an option.  First invoking the functions to compute d1 and d2.  Then returns
the calculation given in the homework description. */
double EuropeanPut(double S, double K, double r, double delta, double sigma, double T, double t)
{
	double _d1 = d1(S, K, r, delta, sigma, T, t);
	double _d2 = d2(_d1, sigma, T, t);
	return ( -S*exp(-delta*(T-t))*F(-_d1) + K*exp(-r*(T-t))*F(-_d2) );
}

/* Computes d1 as given in the Black Scholes equation. */
double d1(double S, double K, double r, double delta, double sigma, double T, double t)
{
	return ((log(S/K)+(r-delta+( pow(sigma,2)/2) )*(T-t))/(sigma*sqrt(T-t)));
}

/* Computes d2 given d1 as a parameter.  The program calculates d1 first, then passes
that value as a parameter to this function. */
double d2(double d1, double sigma, double T, double t)
{
	return (d1 - sigma*(sqrt(T-t)));
}






