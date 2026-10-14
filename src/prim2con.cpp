// Round trip: builds the conserved state of the given primitives, then recovers the
// primitives with both inverters (original: orig.hpp, new: funcs.hpp), to compare with the inputs.
//
// Usage: ./prim2con -rho RHO -P P -lfac W -bb B2 [-cos C] [-tol TOL]
//   -bb   B^2 of the Eulerian field (the same quantity con2prim takes as -bb)
//   -cos  cosine of the angle between v and B (default 0: v perpendicular to B)
//   -tol  solver tolerance (default 0: machine precision)
#include "../include/con2prim/funcs.hpp"
#include "../include/con2prim/orig.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

int main(int argc, char* argv[]) {
    double rho = 0.0;
    double P = 0.0;
    double lfac = 1.0;
    double bb = 0.0;
    double cos_vb = 0.0;
    double tol = 0.0;

    for (int i = 1; i + 1 < argc; i++) {
        if (strcmp(argv[i], "-rho") == 0) {
            rho = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "-P") == 0) {
            P = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "-lfac") == 0) {
            lfac = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "-bb") == 0) {
            bb = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "-cos") == 0) {
            cos_vb = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "-tol") == 0) {
            tol = std::atof(argv[++i]);
        }
    }
    if (!(rho > 0.) || !(P >= 0.) || !(lfac >= 1.) || !(bb >= 0.) || !(std::abs(cos_vb) <= 1.)) {
        fprintf(stderr, "Need rho > 0, P >= 0, lfac >= 1, bb >= 0 and |cos| <= 1\n");
        return 1;
    }

    const double gamma = 5. / 3.;

    // Conserved state, in long double and in cancellation-free forms, rounded once to double.
    // S is split along and across B: S_par = rho h W^2 v_par, S_perp = (rho h W^2 + B^2) v_perp
    using LD = long double;
    const LD g = gamma, r = rho, p = P, W = lfac, B2 = bb, c = cos_vb;
    const LD usq = W * W - 1;             // (W v)^2
    const LD vsq = usq / (W * W);
    const LD v = std::sqrt(vsq);
    const LD u = p / (g - 1);             // internal energy density
    const LD A = (r + u + p) * W * W;     // rho h W^2
    const LD Bm = std::sqrt(B2);
    const LD vB = v * Bm * c;             // v.B
    const LD D = r * W;
    // tau = E - D = rho h W^2 - P + (1 + v^2) B^2 / 2 - (v.B)^2 / 2 - D, without cancellation
    const LD tau = r * W * usq / (W + 1) + u * W * W + p * usq + (1 + vsq) * B2 / 2 - vB * vB / 2;
    const LD Spar = A * v * c;
    const LD Sperp = (A + B2) * v * std::sqrt(1 - c * c);
    const LD SB = A * v * c * Bm;         // S.B

    // Same normalization as con2prim: q = tau/D, r^2 = S^2/D^2, (r.b)^2 = (S.B)^2/D^3, b^2 = B^2/D
    const double d = D;
    const double q = tau / D;
    const double r2 = (Spar * Spar + Sperp * Sperp) / (D * D);
    const double rb2 = SB * SB / (D * D * D);
    const double b2 = B2 / D;

    double P_out, rho_out, lfac_out, xi, mu_max, mu_out;
    auto show = [&](const char* name, const ctp::RootStatus status) {
        std::cout << "--- " << name << " ---"
                  << (status == ctp::RootStatus::Success ? "" : "  (max iterations reached)") << std::endl;
        std::cout << std::setprecision(15);
        std::cout << "P = " << P_out << std::endl;
        std::cout << "rho = " << rho_out << std::endl;
        std::cout << "lfac = " << lfac_out << std::endl;
    };

    show("Original implementation", ctp::get_prims_orig(P_out, rho_out, lfac_out, xi, mu_max, mu_out,
                                        q, d, r2, rb2, b2, gamma, tol, 10000));
    show("\nNew implementation", ctp::get_prims(P_out, rho_out, lfac_out, xi, mu_max, mu_out,
                                   q, d, r2, rb2, b2, gamma, tol, 10000));

    return 0;
}
