/*
Author:		Brian Pendleton
UNCID:		800458057
Date:		11/13/2010
Class:		MATH-6204
Semester:	FALL-2010

Purpose:
This code is designed to numerically solve the heat equation. It will compute the
approximate values of y(x, tao) at x=.5 and tau=.5. 

Coding concepts:
This project uses a class to drive all of the computations. I've made use of typedefs
to allow the code to read much more like math that everything is based on.

Algorithms:
We use the discretized heat equation for the different algorithms noted below:

Explicit Scheme
	Uses a forward difference for the time, and a central difference for the space.
	We simply use a mesh that is [i][j] large, where j=1/.05=20 & i=.5/.001=500

	Builds the mesh using a forward time stepping loop. The logic used is as follows:
	Y[i][j] = lambda*Y[i-1][j-1] + (1-2*lambda)*Y[i-1][j] + lambda*Y[i-1][j+1];


Implicit Scheme (with Thomas Algorithm)
	Uses a backwards difference to calculate the mesh.  We build the mesh using the 
	Thomas algorithm which involves a forward loop that computes b_hat (from slides).
	b_hat is then used in conjunction with the decomposed alpha_hat.

	The RHS of the system of equations is simply:
	b[j] = Y[i-1][j];
		
Crank-Nicolson (with Thomas Algorithm)
	Uses the same algorithm as above, with a slight difference in the math.  This method
	requires the right hand side of the system of equations to be calculated for each time
	step. There are 3 additional terms for the RHS in this method when compared to the 
	Implicit with Thomas Algorithm.

	Here the RHS is calculated according to ...
	b[j] = Y[i-1][j] - gamma[j]*Y[i-1][j-1] + ( gamma[j]+beta[j] )*Y[i-1][j] - beta[j]*Y[i-1][j+1];


Implicit with SOR
	Similar as the Implicit with Thomas, except instead of using a direct solver, we use the 
	iterative technique that stops iterations after a given level of tolerance (error).

Crank-Nicolson with SOR
	Same as Crank-Nicolson with Thomas, except we use the iterative SOR solver instead of the 
	Thomas direct solver.

*/

#if defined(__WIN32__) || defined(_WIN32_) || defined(__WIN32) || defined(_WIN32) || defined(WIN32)
#include "stdafx.h"
#endif
#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include <vector>

using namespace std;

/* Numerical constants */
#define PI 3.14159265358979323846
#define SOR_TOLERANCE 0.00000001

/*
	These are typedefs to make the code read more like the MATH behind the scenes. 
		ARRAY: An array of doubles ( 1 row of our mesh )
		MESH: Our matrix of doubles, implemented as vector of vectors. (array of arrays)
		SUBDIAGONAL: The lower diagonal of the tridiagonal matrix.  (gamma from our slides).
		DIAGONAL: The middle diagonal of the tridiagonal matrix.  (alpha from our slides).
		SUPERDIAGONAL: The upper diagonal of the tridiagonal matrix. (Beta from our slides).
		RHS: The right hand side of our system of equations. (b from our slides)
	
	The typedefs are very useful in the function definitions because 
*/
typedef vector<double> ARRAY;
typedef vector<ARRAY> MESH;
typedef vector<double> SUBDIAGONAL;
typedef vector<double> DIAGONAL;
typedef vector<double> SUPERDIAGONAL;
typedef vector<double> RHS;


/*
	We create a simple class called HeatEquationSolver.
	This class encapsulates all the functionality for this project.  It has solvers for 
	Explicit, Implicit, Crank-Nicolson, and can use techniques such as Thomas and SOR.
*/
class HeatEquationSolver
{
	private:
		double _lambda, _omega, _dx, _dtau;
		bool _sorPrintingEnabled;
		
		// Computational helpers kept private
		ARRAY LUDecomposition(const SUBDIAGONAL &, const DIAGONAL &, const SUPERDIAGONAL &);
		void ThomasSolve(ARRAY &, const SUBDIAGONAL &, const DIAGONAL &, const SUPERDIAGONAL &, const RHS &);
		void SORSolve(ARRAY &, const SUBDIAGONAL &, const DIAGONAL &, const SUPERDIAGONAL &, const RHS &, double);
		void InitializeMesh(MESH &);

		double EuclideanNorm(const ARRAY &, const ARRAY &);
		void PrintSORHeader(string);
		void PrintSORIteration(double, int, double, double);

	public:
		HeatEquationSolver();
		HeatEquationSolver(double, double, double, double);
		
		// Public setters for class properties
		void SetLambda(double);
		void SetOmega(double);
		void SetDeltaX(double);
		void SetDeltaTau(double);
		void EnableSORPrinting();
		void DisableSORPrinting();
		
		// Computational Methods of the class.
		double y(double,double);
		double Explicit(double, double);
		double Implicit_Thomas(double, double);
		double CrankNicolson_Thomas(double, double);
		double Implicit_SOR(double, double);
		double CrankNicolson_SOR(double, double);

		

};

int main(int argc, char* argv[])
{
	cout.sync_with_stdio(false);
	cout << endl;
	HeatEquationSolver eqn = HeatEquationSolver();	// Create a heat equation solver object.
	eqn.EnableSORPrinting();
	//eqn.DisableSORPrinting();  // Uncomment this line to disable SOR Printing (MUCH FASTER EXECUTION TIME). 
	
	double exact = eqn.y(.5, .5);							
	double explicitMethod = eqn.Explicit(.5, .5);			
	double implicitMethod = eqn.Implicit_Thomas(.5, .5);
	double cnMethod = eqn.CrankNicolson_Thomas(.5, .5);
	double cnSORMethod = eqn.CrankNicolson_SOR(.5, .5);		
	double implicitSORMethod = eqn.Implicit_SOR(.5, .5);	

	double exactERR = fabs(exact - exact);
	double explicitERR = fabs(exact - explicitMethod);
	double implicitERR = fabs(exact - implicitMethod);
	double cnERR = fabs(exact - cnMethod);
	double implicitSORERR = fabs(exact - implicitSORMethod);
	double cnSORERR = fabs(exact - cnSORMethod);
	
	// Formatted output for the final results.  The SOR iterations are printed as part of the SOR methods called above.
	cout << "\t\t\t\tSolution Results" << endl << endl;;
	cout << setw(10) << "Method" << setw(18) << "Algorithm" << setw(25) << "Calculated y(.5,.5)" 
		<< setw(35) << "Error ( |Calculated - Exact| )"	<< endl;
	cout << "----------------------------------------------------------------------------------------" << endl ;
	cout.precision(10);
	cout.flags ( ios::fixed );
	cout << setw(18) << left << "EXACT" << setw(15) << "Analytical" << setw(15) << right << exact << setw(30) << exactERR << endl;
	cout << setw(18) << left << "EXPLICIT" << setw(15) << "Forward" << setw(15) << right << explicitMethod << setw(30) << explicitERR << endl;
	cout << setw(18) << left << "IMPLICIT" << setw(15) << "Thomas" << setw(15) << right << implicitMethod << setw(30) << implicitERR << endl;
	cout << setw(18) << left << "CRANK-NICOLOSON" << setw(15) << "Thomas" << setw(15) << right << cnMethod << setw(30) << cnERR << endl;
	cout << setw(18) << left << "IMPLICIT" << setw(15) << "SOR" << setw(15) << right << implicitSORMethod << setw(30) << implicitSORERR << endl;
	cout << setw(18) << left << "CRANK-NICOLOSON" << setw(15) << "SOR" << setw(15) << right << cnSORMethod << setw(30) << cnSORERR << endl;

}

/*
	Default Constructor. Assigns the values according to the project.
*/
HeatEquationSolver::HeatEquationSolver()
{
	_lambda = 0.4; 	_dx = 0.05;
	_dtau = .001;	_omega = 1.2;
	_sorPrintingEnabled = true;
}

/*
	Alternate constructor to assign custom values to the parameters.
*/
HeatEquationSolver::HeatEquationSolver(double lambda, double dx, double dtau, double omega)
{
	_lambda = lambda;	_dx = dx;
	_dtau = dtau;	_omega = omega;
}

/*
	Public setters to modify the properties of the class.
*/
void HeatEquationSolver::SetLambda(double lambda) { _lambda = lambda; }
void HeatEquationSolver::SetOmega(double omega) { _omega = omega; }
void HeatEquationSolver::SetDeltaX(double dx) { _dx = dx; }
void HeatEquationSolver::SetDeltaTau(double dtau) { _dtau = dtau; }
void HeatEquationSolver::EnableSORPrinting() { _sorPrintingEnabled = true; }
void HeatEquationSolver::DisableSORPrinting() { _sorPrintingEnabled = false; }


// Function to calculate the exact value of the heat equation at a given x, tau
double HeatEquationSolver::y(double x, double tau)
{
	return exp( -pow(PI, 2)*tau )*sin( PI*x );
}


double HeatEquationSolver::Explicit(double x, double tau)
{
	int t_steps =  (int)ceil((tau/_dtau)+1);
	int x_steps = (int)ceil((1/_dx)+1);
	int n = x_steps - 1;
	int k = t_steps - 1;
	double lambda2 = (1-(2*_lambda));
	int i, j;

	ARRAY space = ARRAY(x_steps);
	MESH Y = MESH(t_steps, space);		// Creates an [t_steps][x_steps] "mesh" (note typedefs at top)

	InitializeMesh(Y);					// This initializes the mesh with boundary conditions 

	// Populate mesh using the forward difference method
	for ( i=1; i<=k; i++ )
	{
		for ( j=1; j<n; j++ )
		{
			Y[i][j] = _lambda*Y[i-1][j-1] + lambda2*Y[i-1][j] + _lambda*Y[i-1][j+1];
		}
	}
	// This code finds the correct index in the mesh to fetch the value we want to find.
	int x_index = (int)ceil((x/_dx));
	int t_index = (int)ceil((tau/_dtau));
	return Y[t_index][x_index];
	
}


double HeatEquationSolver::Implicit_Thomas(double x, double tau)
{
	int t_steps =  (int)ceil((tau/_dtau)+1);
	int x_steps = (int)ceil((1/_dx)+1);
	int n = x_steps - 1;
	int k = t_steps - 1;
	double lambda2 = (1+(2*_lambda));
	int i;

	ARRAY spaceArray = ARRAY(x_steps);
	MESH Y = MESH(t_steps, spaceArray);
	SUBDIAGONAL gamma = SUBDIAGONAL(x_steps, -_lambda);		
	DIAGONAL alpha = DIAGONAL(x_steps, lambda2);			
	SUPERDIAGONAL beta = SUPERDIAGONAL(x_steps, -_lambda);  

	// Boundary conditions for our coeficients.
	alpha[0] = 1.0;	alpha[n] = 1.0;
	gamma[n] = 0.0;	beta[0] = 0.0;

	InitializeMesh(Y);

	ARRAY alpha_hat = LUDecomposition(gamma, alpha, beta);	// Helper method to create alpha hat.
		
	// Thomas Algorithm
	ARRAY b;
	for ( i=1; i<=k; i++ )
	{
		b = ARRAY(Y[i-1]);
		ThomasSolve(Y[i], gamma, alpha_hat, beta, b);	// Algorithm helper that contains common code for Thomas method. 
	}

	int x_index = (int)ceil((x/_dx));
	int t_index = (int)ceil((tau/_dtau));
	return Y[t_index][x_index];
	
}


double HeatEquationSolver::Implicit_SOR(double x, double tau)
{
	int t_steps =  (int)ceil((tau/_dtau)+1);
	int x_steps = (int)ceil((1/_dx)+1);
	//int n = x_steps - 1;
	int t = t_steps - 1;
	double lambda2 = (1+(2*_lambda));
	

	ARRAY spaceArray = ARRAY(x_steps);
	MESH Y = MESH(t_steps, spaceArray);
	SUBDIAGONAL gamma = SUBDIAGONAL(x_steps, -_lambda);		
	DIAGONAL alpha = DIAGONAL(x_steps, lambda2);			
	SUPERDIAGONAL beta = SUPERDIAGONAL(x_steps, -_lambda);    

	InitializeMesh(Y);

	ARRAY b;
	PrintSORHeader("Implicit SOR Method");
	for ( int i=1; i<=t; i++ )
	{
		b = Y[i-1];
		SORSolve(Y[i], gamma, alpha, beta, b, (i*_dtau));	// SOR Solver that is common code and can be used for Crank-Nicoloson.
	}
	int x_index = (int)ceil((x/_dx));
	int t_index = (int)ceil((tau/_dtau));
	double result = Y[t_index][x_index];
	if ( _sorPrintingEnabled ) { cout << setprecision(10) << fixed << "Final result y(.5, .5) = " << result << endl << endl; }
	return result;

}

// Here we replace _lambda with _lambda/2 because of the equation.
double HeatEquationSolver::CrankNicolson_Thomas(double x, double tau)
{
	int t_steps =  (int)ceil((tau/_dtau)+1);
	int x_steps = (int)ceil((1/_dx)+1);
	int n = x_steps - 1;
	int k = t_steps - 1;
	double lambda2 = 1+_lambda;
	int i, j;

	ARRAY space = ARRAY(x_steps);
	MESH Y = MESH(t_steps, space);
	
	SUBDIAGONAL gamma = ARRAY(x_steps, -.5*_lambda);	// Initialized with -1/2 * lambda
	DIAGONAL alpha = ARRAY(x_steps, lambda2);			// Initialized with 1 + (gamma+beta)
	SUPERDIAGONAL beta = ARRAY(x_steps, -.5*_lambda);   // Initialized with -1/2 * lambda

	// Boundary conditions for our coeficients.
	alpha[0] = 1.0;	alpha[n] = 1.0;
	gamma[n] = 0.0;	beta[0] = 0.0;
	
	InitializeMesh(Y);
	
	ARRAY alpha_hat = LUDecomposition(gamma, alpha, beta);
	
	// Thomas Algorithm
	ARRAY b = ARRAY(x_steps);
	for ( i=1; i<=k; i++ )		
	{
		//b = Y[i-1];
		for ( j=1; j<n; j++ )
		{
			// Crank-Nicolson calculation for the right hand side of the equation
			b[j] = Y[i-1][j] - gamma[j]*Y[i-1][j-1] + ( gamma[j]+beta[j] )*Y[i-1][j] - beta[j]*Y[i-1][j+1];
		}
		ThomasSolve(Y[i], gamma, alpha_hat, beta, b);
	}

	int x_index = (int)ceil((x/_dx));
	int t_index = (int)ceil((tau/_dtau));
	return Y[t_index][x_index];
	
}

double HeatEquationSolver::CrankNicolson_SOR(double x, double tau)
{
	int t_steps =  (int)ceil((tau/_dtau)+1);
	int x_steps = (int)ceil((1/_dx)+1);
	int n = x_steps - 1;
	int k = t_steps - 1;
	double lambda2 = 1+_lambda;
	int i, j;

	ARRAY space = ARRAY(x_steps);
	MESH Y = MESH(t_steps, space);
	
	SUBDIAGONAL gamma = ARRAY(x_steps, -.5*_lambda);	// Initialized with -1/2 * lambda
	DIAGONAL alpha = ARRAY(x_steps, lambda2);			// Initialized with 1 + (gamma+beta)
	SUPERDIAGONAL beta = ARRAY(x_steps, -.5*_lambda);   // Initialized with -1/2 * lambda

	InitializeMesh(Y);
	
	// SOR Algorithm
	ARRAY b = ARRAY(x_steps);
	PrintSORHeader("Crank-Nicolson SOR Method");
	for ( i=1; i<=k; i++ )		
	{
		for ( j=1; j<n; j++ )
		{
			// Crank-Nicolson calculation for the right hand side of the equation
			b[j] = Y[i-1][j] - gamma[j]*Y[i-1][j-1] + ( gamma[j]+beta[j] )*Y[i-1][j] - beta[j]*Y[i-1][j+1];
		}
		SORSolve(Y[i], gamma, alpha, beta, b, (i*_dtau));
	}
	
	
	int x_index = (int)ceil((x/_dx));
	int t_index = (int)ceil((tau/_dtau));
	double result = Y[t_index][x_index];
	if ( _sorPrintingEnabled ) { cout << setprecision(10) << fixed << "Final result y(.5, .5) = " << result << endl << endl; }
	return result;
	
}

// Populates a mesh (passed by reference) with the boundary conditions.
void HeatEquationSolver::InitializeMesh(MESH& Y)
{
	int k = Y.size()-1;
	int n = Y[0].size()-1;

	// Populate mesh with y(x, 0) = sin ( p*x )
	for ( int j=0; j <= n; j++ ) { Y[0][j] = sin( PI*(j*_dx) ); }

	// Boundary conditions for y(0,tau) = y(1,tau) = 0
	for ( int i=0; i <= k; i++ ) { Y[i][0] = 0.0; Y[i][n] = 0.0; }
}


/*
	This function performs an LU decomposition on the coefficients in the tridiagonal matrix.  It creates a 
	new array called alpha_hat (from the slides). This is used in the LUSolve function later when we solve
	the linear systems of equations.
*/
ARRAY HeatEquationSolver::LUDecomposition(const SUBDIAGONAL& gamma, const DIAGONAL& alpha, const SUPERDIAGONAL& beta)
{
	int n = alpha.size()-1;
	ARRAY alpha_hat = ARRAY(n+1);
	alpha_hat[0] = alpha[0];
	for ( int j=1; j<=n; j++ )
	{
		alpha_hat[j] = alpha[j] - beta[j-1]*( gamma[j]/alpha_hat[j-1] );
	}
	return alpha_hat;
}

/*
	This function implements the 2 loops for the Thomas algorith.  It accepts all of the 
	parameters by reference so we avoid needless copies of objects in memeory.  The code
	calculates b_hat (from the slides) using a forward loop, and then goes backwards to find
	the values we're trying to calculate.

	Because our final loop references j+1, at index j, we do not need to create a temporary
	array b_hat.  Instead, we utilize the array x which is our array of X varialbes we're solving
	for.  It was passed by reference so we can just write directly to the mesh.

	The forward loop stores the values in the array, then the backwards loop uses those values 
	and overwrites them going back to index 0.

*/
void HeatEquationSolver::ThomasSolve(ARRAY &x, const SUBDIAGONAL& gamma, const DIAGONAL& alpha_hat, const SUPERDIAGONAL& beta, const RHS &b)
{
	int n = gamma.size()-1;
	int j;

	// Forward loop to calculate b_hat. Values are stored directly in the array x.
	// These values will be used, then overwritten in the backwards loop below.
	x[0] = b[0];
	for ( j=1; j<=n; j++ )
	{
		x[j] = b[j] - x[j-1]*( gamma[j]/alpha_hat[j-1] );
	}

	// Backward loop to calculate Y. 
	x[n] /= alpha_hat[n];		// x_n = b_n / alpha_hat_n
	for ( j=n-1; j>=0; j-- )
	{
		x[j] = (x[j] - (beta[j]*x[j+1]))/alpha_hat[j];		// Calculate x_j and overwrite the values at each step.
	}
	
}

/*
	This method solves the system of equations using the SOR iterative technique.  The algorithm uses a tolerance level
	to stop the iterations of solving for the X vector.  We are able to suppress the internal loop because we're dealing
	with a tridiagonal system. Normaly the algorithm requires you to loop through all of the matrix A when solving, but
	we only need to worry about our 3 vectors.
	And because the subdiagonal is always updated, we know exactly where to use the updated information in our code.

	The solution vector is stored in the array that gets passed by reference.  It is actually mesh[i][] because we
	have fixed our time level and are trying to solve for the space level vector.

*/
void HeatEquationSolver::SORSolve(ARRAY &x, const SUBDIAGONAL& gamma, const DIAGONAL& alpha, const SUPERDIAGONAL& beta, const RHS &b, double tau)
{
	
	int n = gamma.size()-1;
	bool keepIterating = true;
	int k = 1;
	ARRAY oldX;
	while ( keepIterating )
	{
		oldX = x;		
		for ( int j=1; j<n; j++ )
		{
			x[j] = oldX[j] + (_omega/alpha[j])*( b[j] - gamma[j]*x[j-1] - alpha[j]*oldX[j] - beta[j]*oldX[j+1] );
		}
		double norm = EuclideanNorm(x, oldX);		// Calculates the euclidean norm of the two vectors to see if we need to continue.
		if ( tau == .500 ) { PrintSORIteration(tau, k, x[10], norm); }		// Printer helper
		if ( norm <= SOR_TOLERANCE ) { keepIterating = false; }
	
		k++;
	}
}

double HeatEquationSolver::EuclideanNorm(const ARRAY &a, const ARRAY &b)
{
	double norm = 0.0;
	int n = a.size()-1;
	for ( int i=0; i<=n; i++ )
	{
		norm += pow(a[i]-b[i], 2);
	}
	//cout << "    EN = " << sqrt(norm) << endl;
	return sqrt(norm);
}

/* SOR Printing methods */
void HeatEquationSolver::PrintSORHeader(string title)
{
	if ( _sorPrintingEnabled )
	{
		cout << endl;
		cout << title << endl;
		cout << setw(6) << "tau" << setw(6) << "k" << setw(24) << "y(.5, .5)" << setw(30) << "||x(k) - x(k-1)||" << endl;
		cout << "----------------------------------------------------------------" << endl ;
	}
}


void HeatEquationSolver::PrintSORIteration(double tau, int k, double y, double norm)
{
	if ( _sorPrintingEnabled )
	{
		cout << fixed << setprecision(3) << setw(6) << tau << setw(6) << k << setw(25) << right << setprecision(10) << y << setw(25) 
			<< scientific << setprecision(4) << norm << endl;
	}

}




