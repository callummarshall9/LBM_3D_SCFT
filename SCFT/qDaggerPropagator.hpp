//
// Created by callummarshall on 25/11/2019.
//
#include "LBM.hpp"

#ifndef CMDLINE_LBM_QSTARPROPAGATOR_HPP
#define CMDLINE_LBM_QSTARPROPAGATOR_HPP


class qDaggerPropagator : public LBM {
public:
    qDaggerPropagator(int nx, int ny, int nz, std::string velocity_set, double c_s, std::string boundary_condition, double gamma_dot, double N, int N_s, double box_length_rg, field_types field_type, double f, double chiN, double* w_A, double* w_B);
    double field(int x, int y, int z) override;
};


#endif //CMDLINE_LBM_QSTARPROPAGATOR_HPP
