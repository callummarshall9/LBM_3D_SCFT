//
// Created by callummarshall on 25/11/2019.
//

#include "qPropagator.hpp"
#include "qDaggerPropagator.hpp"
#include "../field_types.cpp"
#include <random>
#ifndef CMDLINE_LBM_SCFT_HPP
#define CMDLINE_LBM_SCFT_HPP

class SCFT {
public:
    qPropagator *q_propagator;
    qDaggerPropagator *q_dagger_propagator;
    SCFT(int NX, int NY, int NZ, std::string velocity_set, double c_s, std::string boundary_conditions,
         double gamma_dot, double N, int N_s, double f, double chiN, double box_length_rg,
         field_types field_type, double mixing_parameter, double field_error_threshold, double density_error_threshold, std::string field_initial);
    void Run();
    double Determine_Error();
    void Update_Fields();
    double compute_reduced_density_segment_A(int x, int y, int z, double Q), compute_reduced_density_segment_B(int x, int y, int z, double Q);
    double Determine_Variance_Total();
    void Determine_Density_Differences();
    double compute_Q();
    void perform_timestep();
    void Run_Propagators();
    void output_lbm_data(std::string filename, bool header, bool output_indices);
    void output_max();
    void output_min();
    bool stable();
    bool u_stable();
    bool relaxation_stable();
    void output_parameters();
    void output_field_parameters();
    double find_u_max();
    double tau_LB;


private:
    double delta_x, delta_s;
    double box_length_rg;
    double N;
    int N_s;
    double mixing_parameter,field_error_threshold,variance_threshold;
    double R_g;
    int NX, NY, NZ;
    int s = 0;
    field_types field_type;
    std::string field_initial;
    std::string velocity_set,boundary_conditions;
    double c_s, gamma_dot,f,chiN;
    double *w_A,*w_B,*density_A,*density_B;//FTS fields.
    inline int scalar_index(int x, int y, int z) const {
        return (z * NX * NY) + (y * NX) + x;
    }
    inline int scalar_index(int x, int y, int z, int w) const {
        return (x + y * NX + z * NX * NY + w * NX * NY * NZ);
    }
};


#endif //CMDLINE_LBM_SCFT_HPP
