#ifndef CON2PRIM_TOOLS_HPP
#define CON2PRIM_TOOLS_HPP


#include <cmath>
#include <stdexcept>

namespace ctp{

template<typename T>
inline T abs(T x){
    return x < 0 ? -x : x;
}


enum class Status {
    success,
    max_iter
};


template<typename T, typename Callable>
Status bisect(T& out, Callable&& f, const T& a, const T& b, const T& atol, int max_iter){
    // Uses false position method with Illinois algorithm for acceleration.
    /**
    @param f: a function object which takes a T and returns a T, with a root between a and b
    @param a: lower bound of root
    @param b: upper bound of root
    @param atol: absolute tolerance of f for convergence.
    @param max_iter: maximum number of iterations
    @note If the atol is too small (e.g. 0.), the function will return success if machine precision causes the brackets to collapse before max_iter is reached.

    */
    T zm = a;
    T zp = b;
    T fm = f(zm);
    T fp = f(zp);

    if (fm * fp > 0){
        throw std::runtime_error("Root not bracketed");
    }
    
    // If bracket already within tolerances, return midpoint
    if ((std::abs(zm-zp) < atol) || ((std::abs(fm) + std::abs(fp)) < 2.0*atol)) {
        out = 0.5*(zm + zp);
        return Status::success;
    }

    T z = 0.5*(zm + zp);
    int iter = 0;
    for (; iter < max_iter; iter++) {
        z = (zm*fp - zp*fm)/(fp-fm);  // linear interpolation to point f(z)=0
        T fz = f(z);
        if ( (z-zm) * (z-zp) >= 0) {
            // Check for machine precision (z collapsed to bracket)
            break;
        } else if ((std::abs(zm-zp) < atol) || (std::abs(fz) < atol)) {
            // Check for convergence
            break;
        } else if (fz*fp < 0.0) {
            // assign zm-->zp if root bracketed by [z,zp]
            zm = zp;
            fm = fp;
            zp = z;
            fp = fz;
        } else {
            // assign zp-->z if root bracketed by [zm,z]
            // 1/2 comes from "Illinois algorithm" to accelerate convergence
            fm = 0.5*fm;
            zp = z;
            fp = fz;
        }
    }
    out = std::max(zm, zp); // For our purposes, we want the upper bracket near the solution. max(zm, zp) is required since the illinois algorithm can swap the brackets.
    return iter < max_iter ? Status::success : Status::max_iter;
}





} // namespace ctp



#endif // CON2PRIM_TOOLS_HPP