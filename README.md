#polynomials to interpolate functions

| Feature               | syntax                                                  |
|-----------------------|-------------------------------------------------------- | 
| construction          | `polynomial<> P( {t0,t1,t2,...,tp}, {y0,y1,y2,...,yp} )`|
| evaluation            | `y =  P( time )`                                        |
| differentiation       | `y =  P.deriv( time, order )`                           |
| rootfinding           | `t =  P.findroot( rootfn, tlower, tupper )`             |
| differentiation       | `polynomial<> deriv =  P[order]`                        |
| hermite interpolation | `polynomial<> P( {t0,t0,t0,t1}, {y0,dy0,ddy0,y1} )`     |
| integration           | `P.integrate( t0, t1 )`                                 |
| addition/substraction | `polynomial<> sum = P1 + P2`                            |
| adding points*        | `P.add_pint( tnew, ynew )`                              |

* * limited support for adding derivatives ( that is tnew already exists in t_ )

### plans
* automatic sorting of input, so `t_[0] <= t_[1] <= ... <= t_[p]`
  * in this case `add_point` will have to insert the new point at an appropriate spot
  * this requires to store all finite differences
* add spline class which stores all divided differences for dynamic order interpolation
* set default the identical function for `rootfn` if `state_type` can be compared to `0`. this way, a root
  of the polynomial is found