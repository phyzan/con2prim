#ifndef CON2PRIM_FUNCS_HPP
#define CON2PRIM_FUNCS_HPP


#include "tools.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>

namespace ctp{


template<typename Real>
class KastaunResidual {

    public:

        static constexpr Real GAMMA_MAX = 50.;

        KastaunResidual(Real q, Real d, Real r_sq, Real rb_sq, Real b_sq, Real Gam, Real h0 = 1) : D(d), h0(h0), Gam(Gam), r_sq(r_sq), rb_sq(rb_sq), b_sq(b_sq), q(q) {
            // Cap z0_sq to limit the maximum Lorentz factor (W_max ~ W_max).
            // Without this cap, extreme r_sq (magnetically dominated regions) leads
            // to W >> 1, causing catastrophic cancellation in ehat_mu and garbage
            // primitives that pass the P > 0 / rho > 0 checks without triggering fixup.
            // The cap must be on z0_sq (not v0_sq) so that W_sq/iW_sq, which use
            // 1 + z0_sq as the capped W^2, remain consistent.
            z0_sq   = std::min(r_sq / (h0 * h0), GAMMA_MAX * GAMMA_MAX - 1.0);
            v0_sq   = z0_sq / (1 + z0_sq);
        }

        Real x_mu(const Real mu) const
        {
            return 1 / (1 + mu * b_sq);
        }

        Real rbar_sq(const Real mu, const Real x) const
        {
            return x * (x * r_sq + mu * (1.0 + x) * rb_sq);
        }

        Real qbar_mu(const Real mu, const Real x) const
        {
            const Real mux = mu * x;
            return q - 0.5 * (b_sq + mux * mux * (b_sq * r_sq - rb_sq));
        }

        Real vsq_hat(const Real mu, const Real rbar_sq) const
        {
            Real mu_sq_rbarsq = mu * mu * rbar_sq;
            if (mu_sq_rbarsq < v0_sq){
                return mu_sq_rbarsq;
            } else {
                return v0_sq;
            }
        }

        Real W_sq(const Real mu, const Real rbar_sq) const
        {
            Real vsq_val = mu * mu * rbar_sq;
            if (vsq_val < v0_sq) {
                return 1.0 / (1.0 - vsq_val);
            } else {
                return 1 + z0_sq;
            }
        }

        Real iW_sq(const Real mu, const Real rbar_sq) const
        {
            Real vsq_val = mu * mu * rbar_sq;
            if (vsq_val < v0_sq) {
                return 1.0 - vsq_val;
            } else {
                return 1.0 / (1.0 + z0_sq);
            }
        }

        Real rho_hat(const Real mu, const Real rbar_sq) const
        {
            return D * std::sqrt(iW_sq(mu, rbar_sq));
        }

        Real ehat_mu(const Real mu, const Real rbar_sq_val, const Real q_bar_val) const
        {
            const Real W_val = std::sqrt(W_sq(mu, rbar_sq_val));
            Real kinetic = vsq_hat(mu, rbar_sq_val) * W_sq(mu, rbar_sq_val) / (1.0 + W_val);
            return W_val * (q_bar_val - mu * rbar_sq_val) + kinetic;
        }

        Real Phat(const Real mu) const
        {
            const Real x = x_mu(mu);
            const Real rbar_sq_val = rbar_sq(mu, x);
            const Real qbar_val = qbar_mu(mu, x);
            return D * (Gam - 1.0) * ( (qbar_val - mu * rbar_sq_val) + 1 - std::sqrt(iW_sq(mu, rbar_sq_val)) );
        }

        // Evaluate residual at a value of mu.
        // Kastaun eqn 44
        Real obj_fun(const Real mu) const
        {
            const Real x = x_mu(mu);
            const Real rbarsq = rbar_sq(mu, x);
            const Real qbar = qbar_mu(mu, x);
            const Real iWhat = std::sqrt(iW_sq(mu, rbarsq));
            const Real ehat = std::max(ehat_mu(mu, rbarsq, qbar), 0.);
            // TODO this is ideal-only
            const Real ahat_mod = ehat * (Gam - 1.0);

            // nu_A = h_hat / W = (1 + a_hat)(1 + eps) / W  (paper eq 46)
            const Real nua = (1 + ehat) * (1.0 + ehat + ahat_mod) * iWhat;
            // nu_B = (1 + a_hat) * (1 + qbar - mu*rbarsq)  (paper eq 47)
            const Real nub = (1 + ehat + ahat_mod) * (1.0 + qbar - mu * rbarsq);
            const Real nuhat = std::max(nua, nub); //nua nd nub are multiplied by (1 + eps) here, to avoid the division (...) / (1 + eps) in the paper, which can cause roundoff issues.

            return mu * (nuhat + mu * rbarsq * (1 + ehat)) - (1 + ehat);
        }

        // Residual for finding bracket values
        // Kastaun eqn 49
        Real bound_obj_fun(const Real mu) const
        {
            const Real x = x_mu(mu);
            Real rbar_val = rbar_sq(mu, x);
            return mu * mu * (h0*h0 + rbar_val) - 1;
        }

    // private:
    Real D, h0, Gam; //provided
    Real r_sq, rb_sq, b_sq, q, z0_sq, v0_sq; //derived-normalized
};

template<typename Real>
RootStatus get_prims(Real& P, Real& rho, Real& lfac, Real& xi, Real& mu_max, Real& mu_out, const Real q, const Real D, const Real r_sq, const Real rb_sq, const Real b_sq, const Real Gam, const Real tol, int max_iter) {
    const Real h0 = 1.;
    KastaunResidual res(q, D, r_sq, rb_sq, b_sq, Gam, h0);
    // Find upper bound for mu using eqn 49
    Real mu_min = 0;

    if (res.r_sq < h0 * h0) {
        mu_max = 1.0 / h0;
    } else {
        // Request maximum accuracy for upper bracket (machine precision)
        // as extreme conditions may require that. Usually very few iterations are required anyway,
        // and no EoS calls are needed, so this is not an important bottleneck.
        bisect_illinois<RootEval::Positive>(mu_max, [&res](const Real mu){
            return res.bound_obj_fun(mu);
        }, mu_min, 1.0/h0, 0.0, 0.0, max_iter);
        

        // Per Kastaun et al. Sec IV A: nudge mu_max slightly upward to
        // guarantee that the master function root is strictly contained.
        // For extreme inputs (r^2 >> 1), f(mu_+) is theoretically >= 0
        // but floating-point roundoff can make it slightly negative.
        mu_max = std::min(mu_max * (1.0 + 4.0 * tol) + 4.0 * tol, 1.0 / h0);
    }
    
    Real f_upper = res.obj_fun(mu_max);
    RootStatus result;
    if (f_upper <= 0.) {
        // mu_max is at the root up to machine precision (can happen for extreme
        // inputs where ehat < 0 is clamped).  mu_max ~ root, use directly.
        mu_out = mu_max;
        result = RootStatus::Success;
    } else {
        result = bisect_illinois<RootEval::Positive>(mu_out, [&res](const Real mu){
            return res.obj_fun(mu);
        }, mu_min, mu_max, tol, 0.0, max_iter);
    }

    // Return raw (unclamped) values so the caller can detect unphysical states
    const Real rbar_sq_val = res.rbar_sq(mu_out, res.x_mu(mu_out));
    P = res.Phat(mu_out);
    rho = res.rho_hat(mu_out, rbar_sq_val);
    lfac = std::sqrt(res.W_sq(mu_out, rbar_sq_val));
    // h = res.h_hat(mu_out, rbar_sq_val);
    const Real eps = res.ehat_mu(mu_out, rbar_sq_val, res.qbar_mu(mu_out, res.x_mu(mu_out)));
    const Real ahat = eps * (Gam - 1.0);
    const Real h_hat = 1.0 + eps + ahat;
    xi = h_hat * lfac * D;
    return result;
}


template<typename Real>
void prim2con(Real& d, Real& tau, Real& ss, Real& sb2, const Real P, const Real rho, const Real lfac, const Real BB) {
    Real theta = 1;
    Real gamma = 5. / 3.;
    Real h_times_rho = rho + gamma / (gamma - 1) * P;
    Real z = h_times_rho * lfac * lfac;
    d = rho * lfac;
    Real v2 = 1 - 1 / (lfac * lfac);
    Real SB2 = z*z * v2 * BB * std::cos(theta) * std::cos(theta);
    tau = z - P + BB/2 * (2*z + BB)/(z*z) * SB2;
    ss = 1;
    sb2 = 1;
}


template<typename Real>
void output_result(Real& P, Real& rho, Real& lfac, Real& xi, Real& mu_max, Real& mu_out, const Real q, const Real D, const Real r_sq, const Real rb_sq, const Real b_sq, const Real Gam, const Real tol, int max_iter) {
    RootStatus result = get_prims(P, rho, lfac, xi, mu_max, mu_out, q, D, r_sq, rb_sq, b_sq, Gam, tol, max_iter);
    if (result == RootStatus::Success) {
        std::cout << "Success: ";
    } else if (result == RootStatus::MaxIter) {
        std::cout << "Max iterations reached: ";
    }
    std::cout << std::setprecision(15);
    std::cout << "mu_max = " << mu_max << std::endl;
    std::cout << "mu_out = " << mu_out << std::endl;
    std::cout << "P = " << P << std::endl;
    std::cout << "rho = " << rho << std::endl;
    std::cout << "lfac = " << lfac << std::endl;
    std::cout << "xi = " << xi << " " << (D/mu_out) << std::endl;

}



} // namespace ctp

#endif // CON2PRIM_FUNCS_HPP