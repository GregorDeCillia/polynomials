 
#include<boost/numeric/ublas/vector.hpp>
#include<iostream>
#include<fstream>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

typedef double value_type;					///< define accuracy of values
//typedef value_type time_type;				///< define the accuracy of the time
//typedef vector<value_type> state_type;		///< typically on the bodyvalue_tpye




template <class T>
unbounded_array<T> make_unbounded_array(std::initializer_list<T> list) {
    unbounded_array<T> result(list.size());
    for(unsigned i = 0; i < list.size(); ++ i)
        result[i] = *(list.begin() + i);
    return result;
}








template<class state_type = value_type, class time_type = value_type>
/// polynomials to interpolate functions \f$ f: R\rightarrow R^n\f$
/**
 * 
 * this class gives you
 * 
 * | Feature               | syntax                                                | 
 * |---------------------- | ----------------------------------------------------- | 
 * | construction          | `polynomial P( {t0,t1,t2,...,tp}, {y0,y1,y2,...,yp} )`|
 * | evaluation            | `y =  P( time )`                                      |
 * | differentiation       | `y =  P.deriv( time, oder )`                          |
 * | rootfinding           | `t =  P.findroot( rootfn, tlower, tupper )`           |
 * | differentiation*      | `polynomial deriv =  P[order]`*                       |
 * | hermite interpolation*| `polynomial P( {t0,t0,t0,t1}, {y0,dy0,ddy0,y1} )`     |
 *
 * * * planned 
 * 
 */
class polynomial{
	vector<time_type> t_;							///< the times of evaluation
	vector<state_type> x_;							///< the values of evaluation
	int degree_;									///< the degree
	vector<state_type> newton_coefficients_;		///< newtons coefficients
	vector<value_type> lagrange_coefficients_;		///< lagranges coefficients
	vector<value_type> standard_coefficients_;		///< standard coefficients

	/// Prevent calling the default constructor
	polynomial(){};


	/// calculate the **newton coefficients**. Takes \f$O(p^2)\f$ operations
	/**
	 * calculate the coefficients for the [newton basis](https://en.wikipedia.org/wiki/Newton_polynomial#Definition) i.e. the \f$c_i\f$ in
	 *
	 * \f[ 
	 *     P(t) = \sum_{ i = 0 }^p c_i P_i(t),\ P_0(t) = 1,\ P_i(t) = 
	 \prod_{j=1}^i ( t - t_j )
	 * \f]
	 *
	 * This is used by using the [divided differences formula](https://en.wikipedia.org/wiki/Divided_differences)
	 */
	void calculate_newton_coefficients(){
		newton_coefficients_ = x_;
		for ( int i = 0; i < degree_ + 1; i++ )
			for ( int j = degree_; j > i; j-- )
				newton_coefficients_[j] = ( newton_coefficients_[j] - 
											newton_coefficients_[j-1] )
					/( t_[j] - t_[j-1-i] );
	}

	/// calculate the lagrange coefficients of the polynomias
	void calculate_lagrange_coefficients( time_type t )
	{
		lagrange_coefficients_ = t_;
		for ( int j = 0; j < degree_ + 1; j++ ){
			lagrange_coefficients_[j] = 1.0;
			for ( int m = 0; m < degree_ + 1; m++ ){
				if ( m != j )
					lagrange_coefficients_[j]*= ( t - t_[m] )/( t_[j] - t_[m] );
			}	
		}
	}

	/// evaluate the polynomial using newtons formula
	/**
	 * \f[
	 *      c_0+(t-t_0)\cdot\big(c_1+(t-t_1)\cdot(c_2+...)))\big)
	 * \f]
	 *
	 * Note that this only takes \f$O( p )\f$ oerations
	 */
	state_type evaluate( time_type t )
	{
		state_type x = newton_coefficients_[ degree_ ];
		for ( int i = degree_ - 1; i >= 0; i-- ){
			x *= ( t - t_[ i ] );
			x += newton_coefficients_[ i ];
		}
		return x;
	}

	/// calculates the coefficients a0, a1,a2 for P(t) = a0+a1t+a2t*t+...
	void calculate_standard_coefficients(){
		vector<state_type> derivatives = derivs( 0.0, degree_ );
		standard_coefficients_.resize( degree_ + 1 );
		for ( int i = 0; i < degree_ + 1; i++ ){
			standard_coefficients_[i] = derivatives[i];
			derivatives = derivatives*1/( i + 1 );
		}
	}


public:
	/// Constructor using vectors of evaluation points
	/** 
	 * Suppose, the function to interpolate \f$f\f$ has been 
	 * evaluated at the points
	 *
	 * \f[ (t_0,y_0),(t_1,y_1),...,(t_p,y_p)  \f] 
	 *
	 * i.e. \f$f(t_i) = y_i,\ i = 0,...,p\f$ then this constructs the 
	 * corresponding Polynomial \f$ P\in P_p \f$ with
	 *
	 * \f[ P( t_i ) = y_i,\ i = 0,...,p \f]
	 *
	 * @param[in] t		vector of times where the function has been evaluated
	 *                  \f$ t_0,t_1,...,t_p \f$
	 * @param[in] x		vector of evaluations
	 *                  \f$ x_0,x_1,...,x_p \f$
	 *
	 * the constructor also calculates the coefficients for the newton basis
	 * for later evaluation purposes.
	 */
	polynomial( vector<time_type> t, vector<state_type> x ):
		t_( t ), 
		x_( x ), 
		degree_( t_.size() - 1 )
	{
		if ( t_.size() != x_.size() )
			std::cout << "Error in constructor: dimensions mismatch" << std::endl;
		calculate_newton_coefficients();
	}

	/// Operator, which makes it possible to evaluate the polynomial as \f$P(t)\f$
	/**
	 * Internally, this method uses the **newton basis** for evaluation
	 * and calculates
	 * \f[
	 *      c_0+(t-t_0)\cdot\big(c_1+(t-t_1)\cdot(c_2+...)))\big)
	 * \f]
	 *
	 * Note that this only takes \f$O( p )\f$ oerations
	 */
	state_type operator()( time_type t )
	{
		return evaluate( t );
	}

	/**
	 * This is similar to the last constructor. The number of evaluations
	 * has to be given explicitly as \f$n\f$
	 *
	 */
	/// Alternative constructor using arrays instead of vectors
	polynomial( time_type* t, state_type* x, int n ) 
		: degree_( n - 1 ){
		for ( int i = 0; i < n; i++ ){
			x_[i] = x[i];
			t_[i] = t[i];
		}
		calculate_newton_coefficients();
	}

	/// third constructor of the form P( {t0,t1,t2,..}, {x0,x1,x2,...} )
	polynomial( std::initializer_list<time_type> timelist,
				std::initializer_list<state_type> statelist ):
		x_( make_unbounded_array( statelist ) ),
		t_( make_unbounded_array( timelist ) ),
		degree_( x_.size() - 1 )
	{
		calculate_newton_coefficients();
	}


	/**
	 *
	 * add a new point \f$(t,x)\f$ to the polynomial. This increases
	 * the degree and the *newton coeficcients* get updated. the coefficent
	 * update is of complexicity \f$O(p)\f$
	 *
	 */
	/// add an evaluation to the polynomial. this increases the degree
	void add_point( time_type t, state_type x ){
		if ( x.size() != x_[0].size() )
			std::cout << "Error in add_point: dimensions mismatch" << std::endl;
		// push back t and x
		degree_ ++;
		x_.resize( degree_ + 1 );
		x_( degree_ ) = x;
		t_.resize( degree_ + 1 );
		t_( degree_ ) = t;
		
		// get new coefficient
		newton_coefficients_.resize( degree_ + 1 );
		newton_coefficients_( degree_ ) = x;
		for ( int i = 0; i < degree_ ; i++ )
			newton_coefficients_[ degree_ ] = 
				( newton_coefficients_[ degree_ ] - newton_coefficients_[i] )/
				( t_[ degree_ ] - t_[i] );
	}

	/// alternative evaluation function using the lagrange basis
	/**
	 * This basis is given by
	 * \f[
	 *   P(t) = \sum_{i =0}^p y_iL_i(t),\ L_j( t )= \prod_{m=0,m\neq j}^p\frac{t-t_m}{t_j-t_m}
	 * \f]
	 * The "coefficients" \f$ l_j := L_j(t) \f$ are calculated and then used 
	 * in the above equation. This is numerically not a good method
	 */
	state_type lagrange_evaluate( time_type t )
	{
		calculate_lagrange_coefficients( t );
		state_type x = x_[0]*lagrange_coefficients_[0];
		for ( int i = 1; i < degree_ + 1; i++ )
			x += x_[i]*lagrange_coefficients_[i];
		return x;
	}

	/// calculates the \f$i\f$-th derivtive \f$P^{(i)}(t)\f$ of \f$P\f$ at \f$t\f$ in \f$O(ip)\f$ operations
	/**
	 * the calculation is based on the following recursive formula.
	 * \f{eqnarray*}{
	 *     S_n  (t) &=& c_0     +  (t-t_0)S_{n-1}(t)   \\
	 *     S_n' (t) &=& S_n(t)  +  (t-t_0)S_{n-1}'(t)  \\ 
	 *     S_n''(t) &=& 2S_n(t) +  (t-t_0)S_{n-1}''(t) \\ 
	 *     ...
	 * \f}
	 * Fore more information, read [these lecture notes](http://pages.cs.wisc.edu/~amos/412/lecture-notes/lecture08.pdf)
	 */
	vector<state_type> derivs( time_type t, int i )
	{
		vector<state_type> out( i + 1 );
		vector<state_type> out_bak = out;
		for ( int j = 0; j < i + 1; j++ )
			out[j] = 0*newton_coefficients_[ degree_ ];
		out[0] = newton_coefficients_[ degree_ ];
		for ( int j = degree_-1; j>= 0; j-- ){
			out_bak = out;
			out[0] = newton_coefficients_[ j ] + ( t - t_[ j ] )*out[0];
			for ( int k = 1; k < i+1; k++ ){
				out[k] = k*out_bak[k-1] + ( t - t_[j] )*out[k];
			}
		}
		return out;
	}

	state_type deriv( time_type t, int i ){
		return derivs( t, i )[i];
	}

	/// finds the root of \f$fun( P( t ) )\f$ using the [secant method](https://en.wikipedia.org/wiki/Secant_method)
	/**
	 * the update formula for the time is
	 * \f[ 
	 *     t_{new} = \frac{t_{low} y_{up} - t_{up}y_{low}}{y_{low}-y_{up} }
	 * \f]
	 */
	time_type findroot( value_type(* fun)(state_type) , 
						time_type tlower, time_type tupper )
	{
		value_type ylow = (*fun)( evaluate(tlower) );
		value_type yup  = (*fun)( evaluate(tupper) );
		if ( ylow == 0 ) return tlower;
		else if ( yup == 0 ) return tupper;
		else if ( ylow*yup > 0 ){
			std::cout << "error: no root in the interval" << std::endl;
			return 0;
		}
		if ( ylow > yup ){
			value_type tmp = tlower;
			tlower = tupper;
			tupper = tmp;
		}
		int cnt = 0;
		do{
			time_type tnew = ( tlower*yup - tupper*ylow )/( yup - ylow );
			value_type ynew = (*fun)( evaluate(tnew) );
			if ( ynew > 0 ){
				tupper =  tnew; 
			}
			else{
				tlower = tnew;
			}
		}while( fabs(tupper - tlower) > 1.0e-12 && cnt < 100 );
		return tlower;
	}	

	/// prints the newton formula to the screen
	/**
	 * The result
	 * looks something like
	 * 
	 * `P(t) = 0+(t-0)*(0+(t-1)*(1+(t-17)*0)))`
	 * \f[ P(t) = 0+(t-0)\cdot(0+(t-1)\cdot(1+(t-17)\cdot0))) \f]
	 *
	 */
	void print()
	{
		std::cout << "P(t) = " << newton_coefficients_[0] << "+(t-"<<t_[0]<<")*";
		for ( int i = 0; i < degree_-1; i ++ )
			std::cout << "("<< newton_coefficients_[i+1]<<"+(t-"<<t_[i+1] << ")*";
		std::cout << newton_coefficients_[degree_];
		for ( int i = 0;i < degree_; i++ )
			std::cout << ")";
		std::cout << std::endl;
	}

	polynomial operator[]( int j ){
		// can that be done more efficiently?
		// currently this is an O(p*p) algorithm
		int degree2 = degree_ - j;
		vector<state_type> x2( degree2+1 );
		vector<time_type> t2( degree2+1 );
		for ( int i = 0; i < degree2+1; i++ ){
			t2[i] = t_[i];
			x2[i] = deriv( t_[i], j );
		}
		polynomial P2( t2, x2 );
		return P2;
	}

	state_type standard_evaluate( time_type t ){
		calculate_standard_coefficients();
		state_type sum = standard_coefficients_[ degree_ ];
		for ( int i = 0; i < degree_; i++ ){
			sum *= t;
			sum += standard_coefficients_[ degree_ - 1 - i ];
		}
		return sum;
	}


};
