#include "../include/con2prim/funcs.hpp"
#include <cstdlib>
#include <cstring>
#include <cstdio>


int main(int argc, char* argv[]) {
    double d = 0.0;
    double tau = 0.0;
    double ss = 0.0;
    double sb2 = 0.0;
    double bb = 0.0;
    double tol = 0.0;

    // max_iter is hardcoded for now, but could be made command line argument if needed. h0 is also hardcoded to 1.
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-d") == 0 && i + 1 < argc) {
            d = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "-tau") == 0 && i + 1 < argc) {
            tau = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "-ss") == 0 && i + 1 < argc) {
            ss = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "-sb2") == 0 && i + 1 < argc) {
            sb2 = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "-bb") == 0 && i + 1 < argc) {
            bb = std::atof(argv[++i]);
        } else if (strcmp(argv[i], "-tol") == 0 && i + 1 < argc) {
            tol = std::atof(argv[++i]);
        }
    }

    std::cout << "================================" << std::endl;
    std::cout << "Input parameters:" << std::endl;
    printf("d   = %.15e\n", d);
    printf("tau = %.15e\n", tau);
    printf("ss  = %.15e\n", ss);
    printf("sb2 = %.15e\n", sb2);
    printf("bb  = %.15e\n", bb);
    std::cout << "================================" << std::endl;

    double P, rho, lfac, mu_max, mu_out, xi;


    // rb2, bb, rr, q = SB2/D**3, BB/D, SS/D**2, tau/D
    double rb2 = sb2 / (d * d * d);
    double b2 = bb / d;
    double r2 = ss / (d * d);
    double q = tau / d;

    double gamma = 5. / 3.;


    ctp::output_result(P, rho, lfac, xi, mu_max, mu_out, q, d, r2, rb2, b2, gamma, tol, 100);
    

    return 0;
}
