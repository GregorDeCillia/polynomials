#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE polynomial_class_test

#include <boost/test/unit_test.hpp>
#include <boost/assign/std/vector.hpp>
#include <iostream>

#include "polynomial.hpp"

double tolerance = 1e-10;

BOOST_AUTO_TEST_CASE( test_initialize ) 
{
	// create P(t) = t + 1
	polynomial<> P( {0, 1, 2}, {3, 4, 5 } );

	// test another state_type
	using namespace boost::numeric::ublas;
	typedef float time_type;
	typedef vector<time_type> state_type;
	vector<time_type>  t(5);
	vector< state_type > x(5);
	for ( int i = 0; i < 5; i++ ){
		t[i] = i;
		x[i].resize(2);
		x[i][0] = i;
		x[i][1] = i+2;
	}
	polynomial< state_type, time_type> P2( t, x );
	BOOST_CHECK_CLOSE( P2(10)[1], 10+2, tolerance );
}

// test vector to be used in many subsequent tests
std::vector<double> t {12, 4,13, 8.9, 3.2};

BOOST_AUTO_TEST_CASE( test_evaluate ) 
{
	//create P(t) = t + 3
	polynomial<> P( {0, 1, 2}, {3, 4, 5 } );
	for ( int i = 0; i < t.size(); i++ )
		BOOST_REQUIRE_CLOSE( P( t[i] ), t[i]+3, tolerance );

	// create P2(t) = 4*t*t+6*t+3
	polynomial<> P2( { -1, 0, 1 }, { 1, 3, 13 } );
	// check evaluatios at certain points
	for ( int i = 0; i < t.size(); i++ )
		BOOST_REQUIRE_CLOSE( P2( t[i] ), 4*t[i]*t[i]+6*t[i]+3, tolerance );
}

BOOST_AUTO_TEST_CASE( test_add_point )
{
	// create P(t) = t
	polynomial<> P({0,1},{0,1});
	// extend to P(t) = 2t - t*t
	P.add_point( 2,0 );
	BOOST_REQUIRE_CLOSE( P( 1 + 21 ), P( 1 - 21 ), tolerance );
	BOOST_REQUIRE_CLOSE( P( 1 + 43 ), P( 1 - 43 ), tolerance );

	// give additional condition P'(2)=1
	P.add_point( 0, 1 );
	BOOST_REQUIRE_CLOSE( P.deriv( 0.0, 1 ), 1, tolerance );
	BOOST_REQUIRE_CLOSE( P( 1 ), 1, tolerance );
	BOOST_REQUIRE_CLOSE( P( 2 ), 0, tolerance );
	BOOST_REQUIRE_CLOSE( P( 0 ), 0, tolerance );
}

BOOST_AUTO_TEST_CASE( test_deriv ) 
{
	// create P(t) = t*t
	polynomial<> P( { -1, 0, 1 }, { 1, 0, 1 } );
	// check wether P'(t) = 2*t
	for ( int i = 0; i < t.size(); i++ )
		BOOST_CHECK_CLOSE( P.deriv( t[i], 1 ), 2*t[i], tolerance );
	// check whether P''(t) = 2
	for ( int i = 0; i < t.size(); i++ )
		BOOST_CHECK_CLOSE( P.deriv( t[i], 2 ), 2, tolerance );
	// check whether P'''(t) = 0
	for ( int i = 0; i < t.size(); i++ )
		BOOST_CHECK_CLOSE( P.deriv( t[i], 3 ), 0, tolerance );

	// create P2(t) = 4*t*t+6*t+3
	polynomial<> P2( { -1, 0, 1 }, { 1, 3, 13 } );
	// check wether P2'(t) = 8*t + 6
	for ( int i = 0; i < t.size(); i++ )
		BOOST_CHECK_CLOSE( P2.deriv( t[i], 1 ), 8*t[i] + 6, tolerance );
}

double rootfun( vector< double > x )
{
	return x[0] + x[1];
}

BOOST_AUTO_TEST_CASE( test_find_root )
{
	vector<double> x(3); x[0] = 3; x[1] = 1, x[2] = 13.32;
    vector< vector<double> > x_vec(3); x_vec[0] = x; x_vec[1] = 2*x; x_vec[2] = x;
	vector<double> t(3); t[0] = 0; t[1] = 1; t[2] = 2;
	// create P(t) = ( 2-(t-1)*(t-1) )*[3; 1; 13.32 ]
	// the second * is vector scalar multiplication
	polynomial<vector<double> > P( t, x_vec );
    BOOST_CHECK_CLOSE( P.findroot( rootfun, -10.0, 0.0 ), 1-sqrt(2), tolerance);
    BOOST_CHECK_CLOSE( P.findroot( rootfun, 13.2, 0.2 ), 1+sqrt(2), tolerance);
    BOOST_CHECK_CLOSE( P.findroot( rootfun, -0.2, 12.2 ), 1+sqrt(2), tolerance);
}

BOOST_AUTO_TEST_CASE( test_lagrange_evaluate )
{
	// create polynomial P(t) = ???
	polynomial<> P( {0, 1.3, 2.71}, {3.142, 4.745, 5.125 } );
	// check whether the newton evaluation and the lagrange evaluation give the same result
	for ( int i = 0; i < t.size(); i++ )
		BOOST_CHECK_CLOSE( P( t[i] ), P.lagrange_evaluate( t[i] ), tolerance );
}

BOOST_AUTO_TEST_CASE( test_standard_evaluate )
{
	// create polynomial P(t) = ???
	polynomial<> P( {0, 1.3, 2.71}, {3.142, 4.745, 5.125 } );
	// check whether the newton evaluation and the lagrange evaluation give the same result
	for ( int i = 0; i < t.size(); i++ )
		BOOST_CHECK_CLOSE( P( t[i] ), P.standard_evaluate( t[i] ), tolerance );
}

BOOST_AUTO_TEST_CASE( test_deriv_poly )
{
	// create P(t) = 4*t*t+6*t+3
	polynomial<> P( { -1, 0, 1 }, { 1, 3, 13 } );
	polynomial<> Q = P[1];
	// check whether Q(t) = 8*t + 6
	for ( int i = 0; i < t.size(); i++ )
		BOOST_CHECK_CLOSE( Q( t[i] ), 8*t[i] + 6, tolerance );
	{
		polynomial<> R = P[2];
		// check whether R(t) = 8
		for ( int i = 0; i < t.size(); i++ )
			BOOST_CHECK_CLOSE( R( t[i] ), 8.0, tolerance );
	}
	// check different construction for the second derivative
	{
		polynomial<> R = P[1][0][1];
		for ( int i = 0; i < t.size(); i++ )
			BOOST_CHECK_CLOSE( R( t[i] ), 8.0, tolerance );
	}
}

BOOST_AUTO_TEST_CASE( test_hermite_interpolation )
{
	vector<double> t2( 5 );
	t2[0] = 2.3; t2[1] = 14.53; t2[2] = 13,14; t2[3]= 91.12; t2[4] = -154.1;
	polynomial<> P( 0*t2, t2 );
	for ( int i = 0; i < t2.size(); i++ )
		BOOST_CHECK_CLOSE( P.deriv( 0.0, i ), t2[i], tolerance );
}

BOOST_AUTO_TEST_CASE( test_add )
{
	polynomial<> P1( 3, {-6.1, 2.21, 0.243, 5.3} );
	polynomial<> P2( {2,2.23,-5.23,8.34,-3.1},{7.2,4.32,3.16,-12.1,3.1} );
	polynomial<> P3 = P1 + P2;
	polynomial<> P4 = P1 - P2;
	for ( int i = 0; i < t.size(); i++ ){
		BOOST_REQUIRE_CLOSE( P1(t[i]) + P2(t[i]), P3(t[i]), tolerance );
		BOOST_REQUIRE_CLOSE( P1(t[i]) - P2(t[i]), P4(t[i]), tolerance );
	}
}

BOOST_AUTO_TEST_CASE( integrate )
{
	// create P(t) = (t-1)*(t-1)
	polynomial<> P( {0,1,2}, {1,0,1} );
	BOOST_REQUIRE_CLOSE( P.integrate( 0, 2 ), 2.0/3.0, tolerance );

	// create P(t) = t*t - 2*t + 7
	polynomial<> P2( 0, { 7, -2, 2 } );
	BOOST_REQUIRE_CLOSE( P2.integrate( 0,1 ), 6.0+1.0/3.0, tolerance );
}



