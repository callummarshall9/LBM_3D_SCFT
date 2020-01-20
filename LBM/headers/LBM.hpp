#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include "vector3.hpp"
#ifndef LBM_H_FILE
#define LBM_H_FILE

class LBM {

public:
	LBM(int grid_size, std::string velocity_set, double c_s, double tau, std::string boundary_conditions, double gamma_dot, int runs, double f, double florry_higgs);
	LBM(int nx, int ny, int nz, std::string velocity_set, double c_s, double tau, std::string boundary_condition, double gamma_dot, int runs, double f, double florry_higgs);
	~LBM();
	void set_velocity(int x_field, int y_field, int z_field, double u_x, double u_y, double u_z);//Set velocity at position in velocity field.
	void set_density(int x_field, int y_field, int z_field, double density);//Set density at position in density field.
	double calculate_feq(int i, int j, int k, int w), calculate_feq(int i, int j, int k, int w, double u_le_x);
	void output_density(), output_velocity(), output_indices_file(), output_lbm_data(std::string filename, bool header=true, bool output_indices = false),
            output_f_array(double* f_array, int z_index);//Output functions.
	void set_field_parameters(double N, double Length_Rg, std::string field_type);
	void compute_density_momentum_moment();
	void stream();//Stream the current equilibrium distribution to the next distribution.
	void collision();//Perform the collision step. Assumes delta t / tau = 1.
	void perform_timestep();//Delta t = 1 lattice unit.
	double find_max();//Used to find the max value.
	double find_min();//Used to find the min value.
	void lookup_reverse();//Used to find the negative of the velocity sets.
	virtual double field(int x, int y, int z);
	double compute_reduced_density_segment_plus_variance();
	double* get_chain_propagator();
	int get_time();
	void reset_time();

    int runs;
    double f;
private:
    std::string field_type;
    int time = 0;//Since dt=1, this increments before a 'timestep' calculations are done.
	int NX,NY,NZ;
    double chain_length, box_length_rg, R_g, florry_higgs;
    double c_s,tau_LB;
	vector3<double>* velocity_field;
    vector3<int>* directions;
    int direction_size = 15;
    double *weights, *chainPropagator, *particle_distributions_star,*particle_distributions;//LBM-like quantities.
    inline int scalar_index(int x, int y, int z, int w) const {
        return (x + y * NX + z * NX * NY + w * NX * NY * NZ);
	}
	void output_array(double *array);
	void set_velocity_set(std::string velocity_set);//Used to internally generate velocity_set.
	void initialise();
	//Lattice directions using D3DQ15. assumed speed of sound c_s = 1/sqrt(3).
	int* reverse_indexes;
	std::string boundary_condition;
	double gamma_dot;//Lees-Edwards Shear Rate.
	std::string velocity_set;
	//This will result in a change in the equlibrium function which will be reflected below.
protected:
    inline int scalar_index(int x, int y, int z) const {
return (z * NX * NY) + (y * NX) + x;
}
};

#endif
