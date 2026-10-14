#ifndef CON2PRIM_ORIG_HPP
#define CON2PRIM_ORIG_HPP


#include "tools.hpp"
#include <algorithm>
#include <cstdio>

namespace ctp{

// Floors::LFAC_CAP in KHARMA
static constexpr double LFAC_CAP = 50.;

template<typename Real>
class KastaunResidualOrig
{
  public:
    KastaunResidualOrig(const Real D, const Real q, const Real bsq, const Real bsq_rpsq,
        const Real rsq, const Real rbsq, const Real v0sq, const Real Gam)
        : D_(D), q_(q), bsq_(bsq), bsq_rpsq_(bsq_rpsq), rsq_(rsq), rbsq_(rbsq), v0sq_(v0sq), Gam_(Gam)
    {}

    Real x_mu(const Real mu) const { return 1.0 / (1.0 + mu * bsq_); }
    Real rbarsq_mu(const Real mu, const Real x) const
    {
        return x * (x * rsq_ + mu * (1.0 + x) * rbsq_);
    }
    Real qbar_mu(const Real mu, const Real x) const
    {
        const Real mux = mu * x;
        return q_ - 0.5 * (bsq_ + mux * mux * bsq_rpsq_);
    }
    Real vhatsq_mu(const Real mu, const Real rbarsq) const
    {
        return std::min(mu * mu * rbarsq, v0sq_);
    }
    Real iWhat_mu(const Real vhatsq) const { return std::sqrt(1.0 - vhatsq); }
    Real rhohat_mu(const Real iWhat) const { return D_ * iWhat; }
    Real ehat_mu(const Real mu, const Real qbar, const Real rbarsq, const Real vhatsq,
        const Real What) const
    {
        return What * (qbar - mu * rbarsq) + vhatsq * What * What / (1.0 + What);
    }

    // Kastaun eqn 44
    Real operator()(const Real mu) const
    {
        const Real x = x_mu(mu);
        const Real rbarsq = rbarsq_mu(mu, x);
        const Real qbar = qbar_mu(mu, x);
        const Real vhatsq = vhatsq_mu(mu, rbarsq);
        const Real iWhat = iWhat_mu(vhatsq);
        const Real What = 1.0 / iWhat;
        const Real rhohat = std::max(rhohat_mu(iWhat), 0.);
        const Real ehat = std::max(ehat_mu(mu, qbar, rbarsq, vhatsq, What), 0.);
        const Real Phat = (Gam_ - 1.0) * rhohat * ehat; // ideal gas
        const Real ahat = Phat / (rhohat * (1.0 + ehat));
        const Real nua = (1.0 + ahat) * (1.0 + ehat) * iWhat;
        const Real nub = (1.0 + ahat) * (1.0 + qbar - mu * rbarsq);
        const Real nuhat = std::max(nua, nub);
        const Real muhat = 1.0 / (nuhat + mu * rbarsq);
        return mu - muhat;
    }

    // Kastaun eqn 49
    Real aux_func(const Real mu) const
    {
        const Real x = 1.0 / (1.0 + mu * bsq_);
        const Real rbarsq = x * (rsq_ * x + mu * (1.0 + x) * rbsq_);
        return mu * std::sqrt(1.0 + rbarsq) - 1.0;
    }

  private:
    const Real D_, q_, bsq_, bsq_rpsq_, rsq_, rbsq_, v0sq_, Gam_;
};

// Same signature as get_prims(): P, rho, lfac, xi = h W D, the bound mu_max (= mu_plus) and mu_out
template<typename Real>
RootStatus get_prims_orig(Real& P, Real& rho, Real& lfac, Real& xi, Real& mu_max, Real& mu_out,
    const Real q, const Real D, const Real r_sq, const Real rb_sq, const Real b_sq, const Real Gam,
    const Real tol, int max_iter)
{
    const Real bsq_rpsq = b_sq * r_sq - rb_sq;
    // Cap the solve just past LFAC_CAP (KHARMA then flags W >= LFAC_CAP as bad_lfac)
    const Real v0sq = std::min(r_sq / (1.0 + r_sq), 1.0 - 1.0 / ((LFAC_CAP + 1.) * (LFAC_CAP + 1.)));
    KastaunResidualOrig<Real> res(D, q, b_sq, bsq_rpsq, r_sq, rb_sq, v0sq, Gam);

    mu_max = 1.;
    if (r_sq >= 1.) {
        bisect_illinois<RootEval::Positive>(mu_max, [&res](const Real z) { return res.aux_func(z); },
            (Real)0., (Real)1., tol, tol, max_iter);
    }
    mu_out = mu_max;
    RootStatus result = bisect_illinois<RootEval::Any>(mu_out, [&res](const Real z) { return res(z); },
        (Real)0., mu_max, tol, tol, max_iter);
    if (result == RootStatus::NotBracketed) {
        // Negative at mu_plus only by round-off: the root sits at mu_plus
        mu_out = mu_max;
        result = RootStatus::Success;
    }

    const Real x = res.x_mu(mu_out);
    const Real rbarsq = res.rbarsq_mu(mu_out, x);
    const Real vsq = res.vhatsq_mu(mu_out, rbarsq);
    const Real iW = res.iWhat_mu(vsq);
    const Real W = 1.0 / iW;
    const Real qbar = res.qbar_mu(mu_out, x);
    rho = res.rhohat_mu(iW);
    const Real eps = res.ehat_mu(mu_out, qbar, rbarsq, vsq, W);
    P = (Gam - 1.0) * eps * std::max(rho, 0.); // u = ehat * rho, as in KHARMA
    lfac = W;
    xi = (1.0 + Gam * eps) * W * D; // h W D, h = 1 + eps + P/rho
    return result;
}

} // namespace ctp

#endif // CON2PRIM_ORIG_HPP
