//
// Created by callummarshall on 25/11/2019.
//

#include "qPropagator.hpp"

qPropagator::qPropagator(int NX, int NY, int NZ, std::string velocity_set, double c_s, double tau, std::string boundary_conditions, double gamma_dot, int runs, double f, double florry_higgs,
        double* field_plus, double* field_minus):
    LBM(NX,NY,NZ, velocity_set, c_s, tau, boundary_conditions, gamma_dot, runs, f, florry_higgs), field_plus(field_plus), field_minus(field_minus) {
}

double qPropagator::field(int x, int y, int z) {
    if(get_time() < runs * f) {
        return (-field_plus[scalar_index(x,y,z)] - field_minus[scalar_index(x,y,z)]) / runs;
    } else {
        return (-field_plus[scalar_index(x,y,z)] + field_minus[scalar_index(x,y,z)]) / runs;
    }
    return LBM::field(x,y,z);
}
