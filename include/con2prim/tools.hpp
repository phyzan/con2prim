#ifndef CON2PRIM_TOOLS_HPP
#define CON2PRIM_TOOLS_HPP


#include <cmath>
#include <stdexcept>

namespace ctp{

using std::abs;


enum class RootStatus {
    Success,
    MaxIter,
    NotBracketed
};

enum class RootEval  {
    Positive, // f(root) = 0+
    Negative, // f(root) = 0-
    Any, // No preference
};

template<RootEval RP, typename T, typename Callable>
RootStatus bisect_illinois(
    T&          out,
    Callable&&  f,
    const T     low,
    const T     high,
    const T     ftol,
    const T     xtol,
    const int   max_iter)
{
    // Uses false position method with Illinois algorithm for acceleration.
    /**
    @tparam RP: Side of the root to return: f(out) >= 0, f(out) <= 0, or Any
    @param f: A function takes a T and returns a T, with a root between low and high
    @param low: Lower bound of root
    @param high: Upper bound of root
    @param ftol: Absolute tolerance of f for convergence (0. turns it off)
    @param xtol: Absolute tolerance of the bracket width for convergence (0. turns it off)
    @param max_iter: maximum number of iterations
    @note With ftol and xtol both 0., the function will return success if machine
    precision causes the brackets to collapse before max_iter is reached.

    */
    // The ends swap as the iterations go, so they can be in either order:
    // a is the end kept from earlier iterations, b is the latest iterate
    T a = low;
    T b = high;
    T fa = f(a);
    T fb = f(b);
    T z = 0.5*(a+b);

    // The bracket end on the requested side of the root, or z for Any
    auto side_end = [&]() -> T {
        if constexpr (RP == RootEval::Positive){
            return fa > 0 ? a : b;
        } else if constexpr (RP == RootEval::Negative){
            return fa > 0 ? b : a;
        } else {
            // Round-off can put z just outside
            return std::min(std::max(z, std::min(a, b)), std::max(a, b));
        }
    };

    // Initial checks before entering the main loop
    if (fa * fb > 0){
        return RootStatus::NotBracketed;
    } else if (fa == 0){
        out = a;
        return RootStatus::Success;
    } else if (fb == 0){
        out = b;
        return RootStatus::Success;
    } else if (std::abs(b - a) <= xtol){
        out = side_end();
        return RootStatus::Success;
    }

    // Loop until convergence is met
    // or brackets collapse to machine precision
    for (int iter = 0; iter < max_iter; iter++) {
        z = (a*fb - b*fa)/(fb-fa);  // linear interpolation to point f(z)=0
        if constexpr (RP != RootEval::Any) {
            // Round-off put z on the bracket: one end is at the root, but the end on
            // the requested side may not be. Bisect instead, so that end converges too
            if ((z-a) * (z-b) >= 0) {
                z = 0.5*(a + b);
            }
        }
        if ((z-a) * (z-b) >= 0) {
            // Bracket collapsed to machine precision.
            out = side_end();
            return RootStatus::Success;
        }
        T fz = f(z);

        // Depending on RP, this branch will be compiled away
        // so no need to bloat the code with trivial constexpr branches
        if ((std::abs(fz) <= ftol) &&
                   (RP == RootEval::Any ||
                    (RP == RootEval::Positive ? fz >= 0 : fz <= 0))) {
            // Converged on the requested side of the root
            out = z;
            return RootStatus::Success;
        } else if (fz*fb < 0.0) {
            // Root bracketed by z and b: b becomes the kept end
            a = b;
            fa = fb;
        } else {
            // Root bracketed by a and z: a is kept again
            // Illinois: halve its f to accelerate convergence
            fa = 0.5*fa;
        }
        b = z;
        fb = fz;
        // Converged in x. Tested after the update, so the last iteration still counts
        if (std::abs(b - a) <= xtol){
            out = side_end();
            return RootStatus::Success;
        }
    }

    out = side_end(); // Best guess respecting the requested root position
    return RootStatus::MaxIter;
}


} // namespace ctp



#endif // CON2PRIM_TOOLS_HPP