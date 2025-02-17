#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include "vector3.hpp"
#include "../../field_types.cpp"

#ifndef LBM_H_FILE
#define LBM_H_FILE

class LBM {

public:
    LBM(int nx, int ny, int nz, std::string velocity_set, double c_s, std::string boundary_condition,
            double gamma_dot, double N, int N_s, double box_length_rg, field_types field_type,
            double f, double chiN, double* w_A, double* w_B);
    ~LBM();
    void set_velocity(int x_field, int y_field, int z_field, double u_x, double u_y, double u_z);//Set velocity at position in velocity field.
    void set_density(int x_field, int y_field, int z_field, double density);//Set density at position in density field.
    double calculate_feq(int i, int j, int k, int w), calculate_feq(int i, int j, int k, int w, double u_le_x);
    void output_density(), output_velocity(), output_indices_file(), output_lbm_data(std::string filename, bool header=true, bool output_indices = false),
            output_f_array(double* f_array, int z_index);//Output functions.
    void compute_density_momentum_moment();
    void stream();//Stream the current equilibrium distribution to the next distribution.
    void collision();//Perform the collision step. Assumes delta t / tau = 1.
    void perform_timestep();//Delta t = 1 lattice unit.
    double find_max();//Used to find the max value.
    double find_min();//Used to find the min value.
    void lookup_reverse();//Used to find the negative of the velocity sets.
    virtual double field(int x, int y, int z);
    double* get_chain_propagator();
    void reset_time();
    double find_u_max();
    double f;
    void output_max();
    void output_min();
    void output_field_parameters();
    void output_parameters();
    bool stable();
    bool u_stable();
    bool relaxation_stable();
    double tau_LB;
    double c_s;
private:
    void set_initial_values();
    int viscosity;
    int box_flatten_length;
    int NX,NY,NZ;
    double box_length_rg, R_g, chiN, L,b;
    vector3<double>* velocity_field;
    vector3<int>* directions;
    int direction_size = 15;
    double *weights, *chainPropagator, *particle_distributions_star,*particle_distributions;//LBM-like quantities.
    inline int scalar_index(int x, int y, int z, int w) const {
        return (x + y * NX + z * NX * NY + w * NX * NY * NZ);
    }
    void output_array(double *array);
    void set_velocity_set(std::string velocity_set);//Used to internally generate velocity_set.
    //Lattice directions using D3DQ15. assumed speed of sound c_s = 1/sqrt(3).
    int* reverse_indexes;
    std::string boundary_condition;
    double gamma_dot;//Lees-Edwards Shear Rate.
    std::string velocity_set;
    //This will result in a change in the equlibrium function which will be reflected below.
protected:
    const int delta_x = 1.0;
    double delta_s;
    int s = 0;//Since dt=1, this increments before a 'timestep' calculations are done.
    field_types field_type;
    double *w_A, *w_B;
    inline int scalar_index(int x, int y, int z) const {
        return (z * NX * NY) + (y * NX) + x;
    }
    double N;
    int N_s;
};

#endif
