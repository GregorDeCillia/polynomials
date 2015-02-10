#include "polynomial.hpp"

#include<boost/numeric/ublas/vector.hpp>
#include<iostream>
#include<fstream>                           // create simulation.dat
#include <boost/numeric/ublas/io.hpp>       // print vectors

using namespace boost::numeric::ublas;

typedef double value_type;					///< define accuracy of values
typedef value_type time_type;				///< define the accuracy of the time
typedef vector<value_type> state_type;		///< typicall based on value_type

std::ofstream simFile;

value_type rootfn( state_type x )
{
		return x[1];
}


int main(){
	int degree = 5;
	int dim = 4;

	vector<time_type> t( degree + 1 );
	for ( int i = 0; i < degree + 1; i++ )
		t[i] = i;
	t[2] = 17;

	state_type xo(dim);
	for ( int i = 0; i < dim; i++ )
		xo[i] = i;
	vector<state_type> x( degree + 1 );
	for ( int i = 0; i < degree + 1; i++ )
		x[i] = ( t[i] * t[i] - 1.0 )*xo;

	polynomial<state_type, time_type> P( t, x );
	//	P.add_point( -2.0, 3.0*xo );
	//	P.print();
	for ( double i = 2.5; i < 10; i++ )
		std::cout << i << " " << P(i) << " " 
				  << P.deriv(i, 2) << std::endl;
    
	state_type x2 = P.lagrangeEvaluate( 3.0 );
	std::cout <<  x2[1] << " " << x2[2]  << std::endl;
 
	simFile.open( "simulation.dat" );
	for ( int i = 0; i < 1000; i++ ){
		simFile << i/1000.0 << "\t" << P( i/1000.0 )[1] << std::endl; 
	}
	simFile.close();

	std::cout << "deriv 2: " << P.deriv( 2.0, 1 ) << std::endl;
	std::cout << P( 2.0 ) << std::endl;
	std::cout << P.findroot( rootfn, 0.0, 2.0 ) << std::endl;

	t.resize( 4 );
	vector<time_type> t2 = t;
	for ( int i = 1; i < t.size(); i++ )
		t2[i] = t[i] * t[i] - t[i];
	polynomial<double,double> P2(t,t2);
	std::cout << P2.deriv( 100.0, 1 ) << std::endl;
	P2.print();
}
