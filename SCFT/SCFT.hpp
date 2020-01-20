//
// Created by callummarshall on 25/11/2019.
//

#include "qPropagator.hpp"
#include "qStarPropagator.hpp"

#ifndef CMDLINE_LBM_SCFT_HPP
#define CMDLINE_LBM_SCFT_HPP

class SCFT {
public:
    qPropagator *q_propagator;
    qStarPropagator *q_star_propagator;
    SCFT(int NX, int NY, int NZ, std::string velocity_set, double c_s, double tau, std::string boundary_conditions, double gamma_dot, int runs, double f, double florry_higgs, double chain_length, double length_rg, std::string m_field_type);
    void perform_field_calculations();
    double compute_reduced_density_segment_A(int x, int y, int z), compute_reduced_density_segment_B(int x, int y, int z);
    double compute_reduced_density_segment_plus_variance();
    void update_density_segments();
    double compute_Q();
    void perform_timestep();
    void run_propogators();
    void output_lbm_data(std::string filename, bool header, bool output_indices);
    void output_max();
    void output_min();
    void update_field_minus(double C, double flory_higgs);
    void update_field_plus(double C, double flory_higgs);
private:
    double R_g;
    int NX, NY, NZ;
    int time = 0;
    std::string velocity_set,boundary_conditions;
    double c_s, tau, gamma_dot,f,florry_higgs;
    int runs;
    double *field_plus,*field_minus,*density_A,*density_B;//FTS fields.
    inline int scalar_index(int x, int y, int z) const {
        return (z * NX * NY) + (y * NX) + x;
    }
    inline int scalar_index(int x, int y, int z, int w) const {
        return (x + y * NX + z * NX * NY + w * NX * NY * NZ);
    }
};


#endif //CMDLINE_LBM_SCFT_HPP
