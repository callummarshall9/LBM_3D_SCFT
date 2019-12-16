#include <cmath>
#include <random>
#include <iostream>
#include <filesystem>
#include "LBM.hpp"
#include "SCFT.hpp"

SCFT::SCFT(int NX,int NY, int NZ, std::string velocity_set, double c_s, double tau, std::string boundary_conditions,
           double gamma_dot, int runs, double f, double florry_higgs, double chain_length, double length_rg, std::string field_type) : NX(NX), NY(NY), NZ(NZ), velocity_set(velocity_set),
           c_s(c_s), boundary_conditions(boundary_conditions), gamma_dot(gamma_dot), runs(runs), f(f), florry_higgs(florry_higgs) {
    int box_flatten_length = NX * NY * NZ;
    field_plus = new double[box_flatten_length];
    field_minus = new double[box_flatten_length];
    for(int i = 0; i < box_flatten_length; i++) {
        field_plus[i] = 0.0;
        field_minus[i] = 0.0;
    }
    density_A = new double[box_flatten_length];
    density_B = new double[box_flatten_length];
    q_propagator = new qPropagator(NX, NY, NZ,velocity_set,c_s,tau,boundary_conditions,gamma_dot, runs, f, florry_higgs, field_plus, field_minus);
    q_star_propagator = new qStarPropagator(NX, NY, NZ , velocity_set, c_s, tau, boundary_conditions, gamma_dot, runs, f, florry_higgs, field_plus, field_minus);
    q_propagator->set_field_parameters(chain_length, length_rg, field_type);
    q_star_propagator->set_field_parameters(chain_length, length_rg, field_type);
    const double delta_x = 1.0;//Lattice units.
    const double L = (double)NX * delta_x;
    R_g = L / length_rg;


}

double gauss_noise(double x) {
    std::default_random_engine generator;
    std::normal_distribution<double> dist(0, 1);
    return dist(generator);
}

void SCFT::output_lbm_data(std::string filename, bool header, bool output_indices) {
    std::ofstream output_stream;
    output_stream.open (filename, std::ofstream::out | std::ofstream::app);
    if(header) {
        if(output_indices) {
            output_stream << "x,y,z,";
        }
        output_stream << "p,u_x,u_y,u_z" << '\n';
    }

    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                if(output_indices) {
                    output_stream << x << "," << y << "," << z << ",";
                }
                output_stream << this->q_propagator->get_chain_propagator()[scalar_index(x, y, z, time)] << "," <<
                    this->q_star_propagator->get_chain_propagator()[scalar_index(x,y,z,time)] << "," <<
                    this->density_A[scalar_index(x,y,z)] << "," <<
                    this->density_B[scalar_index(x,y,z)] << "," <<
                    this->q_propagator->field(x,y,z) << "," <<
                    this->q_star_propagator->field(x,y,z) << '\n';
            }
        }
    }
    output_stream.close();

}

void SCFT::output_max() {
    std::cout << "Q propagator max value: " << q_propagator->find_max() << '\n';
    std::cout << "Q dagger propagator max value: " << q_star_propagator->find_max() << '\n';
}

void SCFT::output_min() {
    std::cout << "Q propagator min value: " << q_propagator->find_min() << '\n';
    std::cout << "Q dagger propagator min value: " << q_star_propagator->find_min() << '\n';
}

void SCFT::perform_timestep() {
    time++;
    q_propagator->perform_timestep();
    q_star_propagator->perform_timestep();
}

double SCFT::compute_reduced_density_segment_plus_variance() {
    double sum = 0.0;
    double sumOfSquares = 0.0;
    for(int x = 0; x < NX; x++ ){
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++ ) {
                //double density_segment_A = compute_reduced_density_segment_A(x,y,z);
                //double density_segment_B = compute_reduced_density_segment_B(x,y,z);

                double density_plus = density_A[scalar_index(x,y,z)] + density_B[scalar_index(x,y,z)];
                sum += density_plus;
                sumOfSquares += pow(density_plus, 2.0);
            }
        }
    }
    double mean = sum / (NX * NY * NZ);
    double meanSquares =  sumOfSquares / (NX * NY * NZ);
    return (meanSquares - pow(mean,2.0));
}

double SCFT::compute_Q() {
    double sum = 0.0;
    double* chainPropagator = q_propagator->get_chain_propagator();

    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                sum = sum + chainPropagator[scalar_index(x, y, z, runs)];
            }
        }
    }
    sum /= (NX * NY * NZ);
    return sum;
}


double SCFT::compute_reduced_density_segment_A(int x, int y, int z) {
    double* chainPropagatorQ = q_propagator->get_chain_propagator();
    double* chainPropagatorQStar = q_star_propagator->get_chain_propagator();
    //Integral of q_star * q from 0..f using Simpson 1/3 rule.
    int up_to = f * runs;
    double sum = 0;
    double h = 1;
    for(int s = 0; s < up_to; s++) {
        if(s == 0) {
            sum = sum + h / 3 * chainPropagatorQ[scalar_index(x, y, z, s)] * chainPropagatorQStar[scalar_index(x, y, z, runs - s)];
        } else if(s == up_to) {
            sum = sum + h / 3 * chainPropagatorQStar[scalar_index(x, y, z, s)] * chainPropagatorQStar[scalar_index(x, y, z, runs - s)];
        } else if(s % 2 == 1) {
            sum = sum + h / 3 * 4 * chainPropagatorQ[scalar_index(x, y, z, s)] * chainPropagatorQStar[scalar_index(x, y, z, runs - s)];
        } else if(s % 2 == 0) {
            sum = sum + h / 3 * 2 * chainPropagatorQ[scalar_index(x, y, z, s)] * chainPropagatorQStar[scalar_index(x, y, z, runs - s)];
        }
    }
    sum /= compute_Q();
    return sum;
}

double SCFT::compute_reduced_density_segment_B(int x, int y, int z) {
    double* chainPropagator = q_propagator->get_chain_propagator();
    double* chainPropagatorQStar = q_star_propagator->get_chain_propagator();

    //Integral of q_star * q from fN..N using Simpson 1/3 rule.
    double sum = 0;
    double h = 1;
    int up_to = runs;
    for(int s = f* runs; s < up_to; s++) {
        if(s == 0) {
            sum = sum + h / 3 * chainPropagator[scalar_index(x, y, z, s)] * chainPropagatorQStar[scalar_index(x, y, z, runs - s)];
        } else if(s == up_to) {
            sum = sum + h / 3 * chainPropagator[scalar_index(x, y, z, s)] * chainPropagatorQStar[scalar_index(x, y, z, runs - s)];
        } else if(s % 2 == 1) {
            sum = sum + h / 3 * 4 * chainPropagator[scalar_index(x, y, z, s)] * chainPropagatorQStar[scalar_index(x, y, z, runs - s)];
        } else if(s % 2 == 0) {
            sum = sum + h / 3 * 2 * chainPropagator[scalar_index(x, y, z, s)] * chainPropagatorQStar[scalar_index(x, y, z, runs - s)];
        }
    }
    sum /= compute_Q();
    return sum;
}


void SCFT::update_density_segments() {
    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                density_A[scalar_index(x,y,z)] = compute_reduced_density_segment_A(x,y,z);
                density_B[scalar_index(x,y,z)] = compute_reduced_density_segment_B(x,y,z);
            }
        }
    }
}


void SCFT::run_propogators() {
    time = 0;
    q_propagator->reset_time();
    q_star_propagator->reset_time();
    for(int i = 0; i < runs; i++) {
        perform_timestep();
    }
}

void SCFT::update_field_minus(double C, double flory_higgs) {
    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                double density_minus = density_A[scalar_index(x,y,z)] - density_B[scalar_index(x,y,z)];

                field_minus[scalar_index(x,y,z)] = field_minus[scalar_index(x,y,z)] - C * (-density_minus + 2.0 / flory_higgs * field_minus[scalar_index(x,y,z)]) + gauss_noise(time);
                //field_plus is assumed to always be complex.
            }
        }
    }
}

void SCFT::update_field_plus(double C, double flory_higgs) {
    const double gamma = 1.0;
    int fail_safe = 0;
    const int fail_safe_max = 100;
    while(compute_reduced_density_segment_plus_variance() > 0.0001) {//Ah, crap
        for(int x = 0; x < NX; x++) {
            for (int y = 0; y < NY; y++) {
                for (int z = 0; z < NZ; z++) {
                    double density_plus = density_A[scalar_index(x,y,z)] + density_B[scalar_index(x,y,z)];
                    field_plus[scalar_index(x,y,z)] = field_plus[scalar_index(x,y,z)] - gamma * C * (density_plus - 1.0);
                }
            }
        }
        std::cout << "Relaxing field plus (attempt " << fail_safe + 1 << "/" << fail_safe_max << ")" << '\n';
        run_propogators();
        update_density_segments();
        fail_safe++;
        if(fail_safe > fail_safe_max) {
            std::cout << "ERROR: INFINITE LOOP DETECTED." << '\n';
            exit(-1);
        }
    }
}


void SCFT::perform_field_calculations() {
    const double C = 1.0 * pow(R_g, 3.0) / runs;
    //C=rho_0 * R_g^3/N, Below Equation 9 in text 'Diblock Copolymer Thin Films: A Field-Theoretic Simulation Study
    //by Alexander-Katz and Glenn H. Fredrickson.
    const double flory_higgs = 13.0;

    for(int i = 0; i < runs; i++) {
        std::cout << "Running propagator run (" << (i+1) << "/" << runs << ")" << '\n';
        run_propogators();
        update_density_segments();//Update reduced segment density A and B on each of the lattice nodes.
        update_field_minus(C,flory_higgs);
        update_field_plus(C,flory_higgs);

        output_lbm_data("output/" + std::to_string(i+1)  + ".csv", false, true);

    }

}