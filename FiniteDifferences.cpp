// FiniteDifferences.cpp : Defines the entry point for the console application.
//
#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <cmath>
#include <ctime>
#include <random>
#include <boost/lambda/lambda.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace boost::numeric;
using namespace std;

double BSPDEImplicit(double spot, double rate, double vol, double strike, double barrier, double maturity);
double BSPDECrank_Nicolson(double spot, double rate, double vol, double strike, double barrier, double maturity);
double BSPDEExplicit(double spot, double rate, double vol, double strike, double barrier, double maturity);
double MonteCarloTest(double spot, double rate, double vol, double strike, double barrier, double maturity);

template<class T>
bool InvertMatrix(const ublas::matrix<T>& input, ublas::matrix<T>& inverse)
{
	typedef ublas::permutation_matrix<std::size_t> pmatrix;

	// create a working copy of the input
	ublas::matrix<T> A(input);

	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());

	// perform LU-factorization
	int res = ublas::lu_factorize(A, pm);
	if (res != 0)
		return false;

	// create identity matrix of "inverse"
	inverse.assign(ublas::identity_matrix<T> (A.size1()));

	// backsubstitute to get the inverse
	ublas::lu_substitute(A, pm, inverse);

	return true;
}

int main()
{
	double result = BSPDEImplicit(100,0,0.25,100,130,1);
	//double result = MonteCarloTest(100,0,0.25,100,130,1);
	//double result = BSPDECrank_Nicolson(100,0,0.25,100,130,1);

	cout<<result<<endl;

	system("Pause");

    return 0;
}


///////////////////////// Implicit Finite difference Scheme ///////////////
double BSPDEImplicit(double spot, double rate, double vol, double strike, double barrier, double maturity)
{
	// Grid dimensions
	const int Xn = 32, Tn = 10;

	// Space grid is log spot
	double lnspot = log(spot), lnstrike=log(strike);

	// Grid range is 4 standard deviations apart each side
	double spotMin = lnspot-4*vol*sqrt(maturity), spotMax = lnspot+4*vol*sqrt(maturity);
	double lnbarrier = log(barrier);

	//grid sizes
	double dx = (spotMax-spotMin)/Xn;
	double dt = maturity /Tn;

	ublas::vector<double> spotVec(Xn+1), priceVec(Xn+1);
	ublas::matrix<double> m = ublas::zero_matrix<double>(Xn+1,Xn+1);
	ublas::matrix<double> minv(Xn+1,Xn+1);// = ublas::zero_matrix<double>(Xn+1,Xn+1)

	//specify space grid
	for(int i=0;i<Xn+1;++i)
	{
		spotVec[i]=spotMin+i*dx;
	}

	//Up and Out Barrier Call payoff at expiry
	for(int i=0;i<Xn+1;++i)
	{
		if(spotVec[i]<lnbarrier)
			priceVec[i]=max(exp(spotVec[i])-exp(lnstrike),0.0);
		else
			priceVec[i]=0.0;
	}

	// Coefficients of the linear equations
	double c1,c2,c3;
	c1 = (rate-0.5*pow(vol,2))*dt/(2*dx)-0.5*pow(vol,2)*dt/pow(dx,2);
	c2 = 1+rate*dt+pow(vol,2)*dt/pow(dx,2);
	c3 = -((rate-0.5*pow(vol,2))*dt/(2*dx)+0.5*pow(vol,2)*dt/pow(dx,2));

	// Create tri-diagonal matrix which will be state independent
	m(0,0) = 1+rate*dt/dx+rate*dt;
	m(0,1) = -rate*dt/dx;
	m(Xn,Xn-1) = rate*dt/dx;
	m(Xn,Xn) = 1-rate*dt/dx+rate*dt;

	for(int a=1;a<Xn;++a)
	{
		m(a,a-1)=c1;
		m(a,a)=c2;
		m(a,a+1)=c3;
	}

	//Invert the matrix
	InvertMatrix(m,minv);

	//Backward induction to arrive at price at time zero
	for (int j=0;j<Tn;++j)
	{
		priceVec = ublas::prod(minv,priceVec);
		for(int i=0;i<Xn+1;++i)
		{
			if(spotVec[i]>lnbarrier)
				priceVec[i]=0.0;
		}
		cout<<priceVec[Xn/2]<<endl;
	}

	//Desired option price
	return priceVec[Xn/2];
}


///////////////////////// Crank-Nicolson Finite difference Scheme ///////////////
double BSPDECrank_Nicolson(double spot, double rate, double vol, double strike, double barrier, double maturity)
{
		// Grid dimensions
	const int Xn = 32, Tn = 10;

	// Space grid is log spot
	double lnspot = log(spot), lnstrike=log(strike);

	// Grid range is 4 standard deviations apart each side
	double spotMin = lnspot-4*vol*sqrt(maturity), spotMax = lnspot+4*vol*sqrt(maturity);
	double lnbarrier = log(barrier);

	//grid sizes
	double dx = (spotMax-spotMin)/Xn;
	double dt = maturity /Tn;

	ublas::vector<double> spotVec(Xn+1), priceVec(Xn+1), tempVec(Xn+1);
	ublas::matrix<double> m = ublas::zero_matrix<double>(Xn+1,Xn+1);
	ublas::matrix<double> tempMatrix(Xn+1,Xn+1);
	ublas::matrix<double> q = ublas::zero_matrix<double>(Xn+1,Xn+1);
	ublas::matrix<double> qinv(Xn+1,Xn+1);

	//specify space grid
	for(int i=0;i<Xn+1;++i)
	{
		spotVec[i]=spotMin+i*dx;
	}

	//Up and Out Barrier Call payoff at expiry
	for(int i=0;i<Xn+1;++i)
	{
		//if(spotVec[i]<lnbarrier)
			priceVec[i]=max(exp(spotVec[i])-exp(lnstrike),0.0);
		//else
		//	priceVec[i]=0.0;
	}

	// Coefficients of the linear equations
	double c1,c2,c3;
	c1 = -0.25*(rate-0.5*pow(vol,2))*dt/dx+0.25*pow(vol,2)*dt/pow(dx,2);
	c2 = -0.5*rate*dt-0.5*pow(vol,2)*dt/pow(dx,2);
	c3 = 0.25*(rate-0.5*pow(vol,2))*dt/dx+0.25*pow(vol,2)*dt/pow(dx,2);

	// Create tri-diagonal matrices which will be state independent
	m(0,0) = 1-0.5*(rate-0.5*pow(vol,2))*dt/dx-0.5*rate*dt;
	m(0,1) = 0.5*(rate-0.5*pow(vol,2))*dt/dx;
	m(Xn,Xn-1) = -0.5*(rate-0.5*pow(vol,2))*dt/dx;
	m(Xn,Xn) = 1+0.5*(rate-0.5*pow(vol,2))*dt/dx-0.5*rate*dt;

	for(int a=1;a<Xn;++a)
	{
		m(a,a-1)=c1;
		m(a,a)=1+c2;
		m(a,a+1)=c3;
	}

	// Now matrix q
	q(0,0) = 1+0.5*(rate-0.5*pow(vol,2))*dt/dx+0.5*rate*dt;
	q(0,1) = -0.5*(rate-0.5*pow(vol,2))*dt/dx;
	q(Xn,Xn-1) = 0.5*(rate-0.5*pow(vol,2))*dt/dx;
	q(Xn,Xn) = 1-0.5*(rate-0.5*pow(vol,2))*dt/dx+0.5*rate*dt;

	for(int a=1;a<Xn;++a)
	{
		q(a,a-1)=-c1;
		q(a,a)=1-c2;
		q(a,a+1)=-c3;
	}


	// Invert matrix q
	InvertMatrix(q,qinv);

	// Multiplying out square matrices
	tempMatrix = ublas::prod(q,m);

	//Backward induction to arrive at price at time zero
	for (int j=0;j<Tn;++j)
	{
		priceVec = ublas::prod(tempMatrix,priceVec);
		/*for(int i=0;i<Xn+1;++i)
		{
			if(spotVec[i]>lnbarrier)
				priceVec[i]=0.0;
		}*/
		cout<<priceVec[Xn/2]<<endl;
	}

	//Desired option price
	return priceVec[Xn/2];
}


///////////////////////// Explicit Finite difference Scheme ///////////////
double BSPDEExplicit(double spot, double rate, double vol, double strike, double barrier, double maturity)
{
	//Grid demensions
	int Xn=20, Tn=10;

	// Space grid is logspot
	double lnspot = log(spot), lnstrike = log(strike);

	// Grid range is 4 standard deviations apart each side
	double spotMin = lnspot-4*vol*sqrt(maturity), spotMax = lnspot+4*vol*sqrt(maturity);
	double lnbarrier = log(barrier);

	//Grid sizes
	double dx = (spotMax-spotMin)/Xn;
	double dt = maturity/Tn;

	//Vectors for space, price, the first derivative & second derivative of price
	vector<double> spotVec,priceVec(Xn+1),diff1Vec(Xn+1),diff2Vec(Xn+1);

	//Specify space grid
	for(int i=0;i<Xn+1;++i)
	{
		spotVec.push_back(spotMin+i*dx);
	}

	//Call option payoff at expiry
	for(int i=0;i<Xn+1;++i)
	{
		if(spotVec[i]<lnbarrier)
			priceVec[i]=max(exp(spotVec[i])-exp(lnstrike),0.0);
		else
			priceVec[i]=0.0;
	}

	//Backward induction from expiry to start
	for(int j=0;j<Tn;++j)
	{
		//First derivative vector
		diff1Vec[0]=(priceVec[1]-priceVec[0])/dx;
		for(int i=0;i<Xn;++i)
		{
			diff1Vec[i]=(priceVec[i+1]-priceVec[i-1])/(2*dx);
		}
		diff1Vec[Xn]=(priceVec[Xn]-priceVec[Xn-1])/dx;

		//Second derivative vector
		diff2Vec[0]=0;
		for(int i=0;i<Xn;++i)
		{
			diff2Vec[i]=(priceVec[i+1]+priceVec[i-1]-2*priceVec[i])/dx*dx;
		}
		diff2Vec[0]=0;

		//Fill in barrier price vector
		for(int i=0;i<Xn+1;++i)
		{
			if(spotVec[i]>lnbarrier)
				priceVec[i]=0.0;
			else
				priceVec[i] = (priceVec[i]+(rate-0.5*vol*vol)*diff1Vec[i]*dt+0.5*pow(vol,2)*diff2Vec[i]*dt)/(1+rate*dt);
		}
	}

	//Desired barrier option price
	return priceVec[Xn/2];
}


///////////////////////// Monte-carlo simulation to test FDM schemes ///////////////
double MonteCarloTest(double spot, double rate, double vol, double strike, double barrier, double maturity)
{
	int numSim = 10000, steps = 50;
	double simSpot = spot, dt = maturity/steps;
	double payoff, sumpayoff=0.0;
	int barrHit = 0;

	// Generate random numbers
	std::default_random_engine rnd(time(0));
	normal_distribution<double> distribution(0,1);
	
	// Simulate paths and get payoffs
	for(int i=0;i<numSim;++i)
	{
		for(int j=0;j<steps;++j)
		{
			simSpot = simSpot*exp((rate-0.5*vol*vol)*dt+vol*sqrt(dt)*distribution(rnd));
			if(simSpot > barrier)
				barrHit = 1;
		}

		// Up and Out barrier call payout
		payoff = (1-barrHit)*max(simSpot-strike,0.0);
		sumpayoff = sumpayoff + payoff/numSim;
		barrHit = 0, simSpot = spot;
	}
	
	return sumpayoff*exp(-rate*maturity);
}












