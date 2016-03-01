/*
Author: Brian Pendleton
ID: 800458057

To build the tree, I used two classes.  One class represents the nodes of the tree,
the other class represents the J&R tree itself.

I used a VECTOR of VECTORS of tree nodes to represent my mesh.

*/

#if defined(__WIN32__) || defined(_WIN32_) || defined(__WIN32) || defined(_WIN32) || defined(WIN32)
#include "stdafx.h"
#endif
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

using namespace std;

#if !defined(MATH_CONSTANTS)
#define MATH_CONSTANTS
#define ERR 4.430465
#endif

class TreeNode
{
	private:
		double _stockPrice;
		double _optionPrice;

	public:
		TreeNode();
		double S();
		void S(double);
		double V();
		void V(double);
		double Error();
};


class JRBinomialTree
{
	private:
		double _initStockPrice;
		double _strike;
		double _T;
		double _r;
		double _sigma;
		double _dt;
		double _beta;
		double _u;
		double _d;
		double _q;
		int _M;
		vector<vector<TreeNode> > _tree;  // This is the mesh. It is a matrix of tree nodes.

		void CalculateDeltaT();
		void CalculateBeta();
		void CalculateU();
		void CalculateD();
		void CalculateQ();
		void CalculateParameters();
		TreeNode GetNode(int, int);

	public:
		JRBinomialTree(double, double, double, double, double);
		void CalculateEuropeanPut(int);
		void PrintCalculatedParameters();
		void PrintResults();
		void BuildTree();
		double CalculatePrice(int, int);
		double PutPayoff(double);
		double RiskNeutralOptionValue(double, double);
		void EuropeanPut();

};




int main(int argc, char* argv[])
{
	double r = 0.06;
	double sigma = 0.3;
	double T = 1.0;
	double K = 10;
	double S = 5;
	
	cout << "Number of Steps (M)" << setw(28) << "Initial Price (V)" << setw(28) << "Error (V - 4.430465)" << endl;
	cout << "------------------------------------------------------------------------------------" << endl;
			
	JRBinomialTree tree = JRBinomialTree(S, K, T, sigma, r);
	tree.CalculateEuropeanPut(5);
	tree.CalculateEuropeanPut(10);
	tree.CalculateEuropeanPut(20);
	tree.CalculateEuropeanPut(50);
	tree.CalculateEuropeanPut(100);
	tree.CalculateEuropeanPut(200);
	tree.CalculateEuropeanPut(500);
	tree.CalculateEuropeanPut(1000);

}



TreeNode::TreeNode() 
{
	_stockPrice = 0.0;
	_optionPrice = 0.0;
}

double TreeNode::S() { return _stockPrice; }
void TreeNode::S(double stockPrice) { _stockPrice = stockPrice; }
double TreeNode::V() { return _optionPrice; }
void TreeNode::V(double optionPrice) { _optionPrice = optionPrice; }
double TreeNode::Error() { return _optionPrice - ERR; }



void JRBinomialTree::CalculateDeltaT() { _dt = _T/_M; }

void JRBinomialTree::CalculateBeta() { _beta = .5*(exp(-_r*_dt)+exp((_r+pow(_sigma,2))*_dt)); }

void JRBinomialTree::CalculateU() { _u = _beta + sqrt(pow(_beta,2)-1); }

void JRBinomialTree::CalculateD() { _d = _beta - sqrt(pow(_beta,2)-1); }

void JRBinomialTree::CalculateQ() { _q = (exp(_r*_dt)-_d)/(_u-_d); }

void JRBinomialTree::CalculateParameters()
{
	CalculateDeltaT();
	CalculateBeta();
	CalculateU();
	CalculateD();
	CalculateQ();
}

TreeNode JRBinomialTree::GetNode(int i, int j) { return _tree[i][j]; }

JRBinomialTree::JRBinomialTree(double S, double K, double T, double sigma, double r)
{
	_initStockPrice = S;
	_strike = K;
	_T = T;
	_r = r;
	_sigma = sigma;
}

void JRBinomialTree::CalculateEuropeanPut(int M)
{
	_M = M;
	_tree.clear();
	CalculateParameters();
	BuildTree();
	EuropeanPut();
	PrintResults();
}

void JRBinomialTree::PrintResults()
{
	TreeNode initialNode = GetNode(0,0);
	cout.precision(6);
	cout << setw(8) << _M << setw(30) << initialNode.V() << setw(34) << setprecision(6) <<  initialNode.Error() << endl;
}

void JRBinomialTree::BuildTree()
{
	double stockPrice;
	
	_tree = vector<vector<TreeNode> >(_M+1);
	for ( int i=0; i<=_M; i++ )
	{
		vector<TreeNode> nodes = vector<TreeNode>(i+1);
		for ( int j=0; j<=i; j++ )
		{
			TreeNode tn = TreeNode();
			stockPrice = CalculatePrice(i,j);
			tn.S(stockPrice);
			nodes[j] = tn;
		}
		_tree[i] = nodes;
	}

}

double JRBinomialTree::CalculatePrice(int i, int j)
{
	int upMoves = j;
	int downMoves = i-j;
	return _initStockPrice*(pow(_u, upMoves)*pow(_d,downMoves));
}

double JRBinomialTree::PutPayoff(double S)
{
	return max(0.0, (_strike-S));
}

double JRBinomialTree::RiskNeutralOptionValue(double vUp, double vDown)
{
	return exp(-_r*_dt)*(_q*vUp + (1-_q)*vDown);
}

void JRBinomialTree::EuropeanPut()
{
	double optionPrice;
	for ( int i=_M; i>=0; i-- )
	{
		for ( int j=0; j<=i; j++ )
		{
			if ( i == _M )
			{
				optionPrice = PutPayoff(_tree[i][j].S());
			}
			else
			{
				TreeNode upMove = _tree[i+1][j+1];
				TreeNode downMove = _tree[i+1][j];
				optionPrice = RiskNeutralOptionValue(upMove.V(), downMove.V());
			}
			_tree[i][j].V((optionPrice));
		}
	}
}

