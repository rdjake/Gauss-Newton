#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <random>
#include <cmath>
#include <stdio.h>
#include <algorithm>
#include <direct.h>
#include <fstream>
#include <string>
#include <cmath>
#include <omp.h>
#include <filesystem>
#include "windows.h"

using namespace std;
typedef vector<double> Vector;
typedef vector<Vector> Matrix;

// К мД = K м2 * 1Е15  | 400мД = 4Е-13
// K = (10мД , 5000 мД), Phi = (0,1) 

class GaussNewton {
	
public :
	int ITERATIONS, n, m, N, M, numOfWells, numOfLayers, numOfContours, timePoints, Threads, numOfPlots;
	string path, pathMesh, pathOutput, pathReverse, pathOutput3D, pathHome, pathProperties, trash;
	Vector TrueParams, Params, TrueFx, Fx;

	vector<bool> MaskForLayers, MaskForWells; //маски должны быть полноразмерные
	vector<bool>  WellsType; //true - добывающие / false - нагнетающие
	vector<vector<string>> plots;
	vector<vector<int>> ContourBound;
	bool PARALLEL, newWay, seeOutput, seeJacobian, useInWells, useOutWells, useS, useM, useP, copyAll, customFx;
	int MaxIters = 1000, numOfIn = 0, numOfOut = 0, offLayerSize, offContoursSize;
	double Func, Func0, Eps = 1e-13, eps = 1e-30, BetaMin = 1e-4, AlphaMin = 1e-15, AlphaMax = 1e+15;
	Vector p_init, p0, p, Delta_k, f_, alpha, falpha;
	Matrix A, Aalpha, J, FunctionalProgress;
	double  scale = 1E13, IntervalK[2] = { 10.0 * scale / 1E15, 5000.0 * scale / 1E15 }, IntervalPhi[2] = { 1E-2, 1.0 };

	inline int sign(double val);

	void createPath();
	void ReadData();
	void copyFiles(string from, string to);
	void SetReverse();
	void readFx(Vector& fx, string pp);
	void setParameters(Vector& params, string pp);
	void init(double k, double phi);
	
	inline void f(Vector& params, Vector& fx);
	double Functional(Vector& truefx, Vector& fx);
	void FunctionalSeparate(Vector& truefx, Vector& fx, int iter);

	void calcJ();
	Vector SolveGauss(Matrix& A, Vector& b);	
	void method();
	
	void Results();
	void ReturnTrue(Vector& params, string pp);
}; 
