//
// Created by callummarshall on 25/11/2019.
//

#include "LBM.hpp"

#ifndef CMDLINE_LBM_QPROPAGATOR_HPP
#define CMDLINE_LBM_QPROPAGATOR_HPP


class qPropagator : public LBM {
public:
    qPropagator(int NX, int NY, int NZ, std::string velocity_set, double c_s, double tau, std::string boundary_conditions, double gamma_dot, int runs, double f, double florry_higgs,
                double* field_plus, double* field_minus);
    double field(int x, int y, int z) override;
    double *field_plus,*field_minus;
};


#endif //CMDLINE_LBM_QPROPAGATOR_HPP
