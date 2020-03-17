#include <cmath>
#include <random>
#include <iostream>
#include <algorithm>
#include <filesystem>
#include "LBM.hpp"
#include "SCFT.hpp"

SCFT::SCFT(int NX, int NY, int NZ, std::string velocity_set, double c_s, std::string boundary_conditions,
           double gamma_dot, double N, int N_s, double f, double chiN, double box_length_rg,
           std::string m_field_type, double mixing_parameter) : NX(NX), NY(NY), NZ(NZ), velocity_set(velocity_set),
                                                                          c_s(c_s), boundary_conditions(boundary_conditions), gamma_dot(gamma_dot), N(N), N_s(N_s), f(f), box_length_rg(box_length_rg), chiN(chiN), field_type(m_field_type), mixing_parameter(mixing_parameter) {
    int box_flatten_length = NX * NY * NZ;
    w_A = new double[box_flatten_length];
    w_B = new double[box_flatten_length];
    this->delta_x = 1.0;//Lattice units.
    this->N = N;
    this->box_length_rg = box_length_rg;
    const double L = (double)NX * delta_x;
    this->R_g = L / box_length_rg;
    const double b = sqrt(6.0) / sqrt(N) * R_g;
    this->delta_s = 1.0;
    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                const double value = 1.0 * sin(2 * M_PI * (x * delta_x - L / 2.0) / L);
                w_A[scalar_index(x,y,z)] = -value;
                w_B[scalar_index(x,y,z)] = value;
            }
        }
    }

    /*

    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                //w_A[scalar_index(x,y,z)] = dist(generator) /N * 10.0;
                w_A[scalar_index(x,y,z)] = 0.0;
                w_B[scalar_index(x,y,z)] = dist(generator) / N;
                //w_A[scalar_index(x,y,z)] = 0.0;
                //w_B[scalar_index(x,y,z)] = 0.0;
            }
        }
    }*/

    //w_B[scalar_index(0,0,0)] = -5.0;
    //w_B[scalar_index(midpoint, 0, 0)] = 5.0;

    density_A = new double[box_flatten_length];
    density_B = new double[box_flatten_length];
    q_propagator = new qPropagator(NX,NY,NZ,velocity_set,c_s,boundary_conditions,gamma_dot,N,N_s,box_length_rg,field_type,f,chiN,w_A,w_B);
    q_dagger_propagator = new qDaggerPropagator(NX, NY, NZ, velocity_set, c_s, boundary_conditions, gamma_dot, N, N_s, box_length_rg, field_type, f, chiN, w_A, w_B);
    tau_LB = q_propagator->tau_LB;
    std::cout << "chiN parameter: " << chiN << '\n';
}



void SCFT::output_lbm_data(std::string filename, bool header, bool output_indices) {
    std::ofstream output_stream;
    output_stream.open (filename, std::ofstream::out | std::ofstream::app);
    if(header) {
        if(output_indices) {
            output_stream << "x,y,z,";
        }
        output_stream << "q,q*,phi_A,phi_B,phi, wA,wB" << '\n';
    }
    std::cout << "s: " << s << '\n';
    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                if(output_indices) {
                    output_stream << x << "," << y << "," << z << ",";
                }
                double phi = density_A[scalar_index(x,y,z)] + density_B[scalar_index(x,y,z)];
                output_stream << this->q_propagator->get_chain_propagator()[scalar_index(x, y, z, s)] << "," <<
                              this->q_dagger_propagator->get_chain_propagator()[scalar_index(x, y, z, s)] << "," <<
                              this->density_A[scalar_index(x,y,z)] << "," <<
                    this->density_B[scalar_index(x,y,z)] << "," <<
                    phi << "," <<
                    w_A[scalar_index(x,y,z)] << "," <<
                    w_B[scalar_index(x,y,z)] << '\n';
            }
        }
    }
    output_stream.close();

}

void SCFT::output_max() {
    std::cout << "Q propagator max value: " << q_propagator->find_max() << '\n';
    std::cout << "Q dagger propagator max value: " << q_dagger_propagator->find_max() << '\n';
}

void SCFT::output_min() {
    std::cout << "Q propagator min value: " << q_propagator->find_min() << '\n';
    std::cout << "Q dagger propagator min value: " << q_dagger_propagator->find_min() << '\n';
}

void SCFT::perform_timestep() {
    s++;
    q_propagator->perform_timestep();
    q_dagger_propagator->perform_timestep();
}

double SCFT::Determine_Variance_Total() {
    double sum = 0.0;
    double sumOfSquares = 0.0;
    for(int x = 0; x < NX; x++ ){
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++ ) {
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
    double sum_q = 0.0;
    double sum_q_star = 0.0;
    double* chainPropagator = q_propagator->get_chain_propagator();
    double* qStarChainPropagator = q_dagger_propagator->get_chain_propagator();
    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                sum_q = sum_q + chainPropagator[scalar_index(x, y, z, N)];
                sum_q_star = sum_q_star + qStarChainPropagator[scalar_index(x,y,z,N)];
            }
        }
    }
    sum_q /= (NX * NY * NZ);
    sum_q_star /= (NX * NY * NZ);

    if((sum_q - sum_q_star) > 0.0001) {//If the difference in the q and q star Q value propagations don't match, output a warning.
        //std::cout << "Warning: Q for q propagator doesn't match q star propagator." << '\n';
    }
    if(sum_q > sum_q_star) {
        return sum_q;
    } else {
        return sum_q_star;
    }
}


double SCFT::compute_reduced_density_segment_A(int x, int y, int z) {
    double* chainPropagatorQ = q_propagator->get_chain_propagator();
    double* chainPropagatorQStar = q_dagger_propagator->get_chain_propagator();
    double sum = 0.0;
    const int up_to = lround(f * N);
    for(int s = 0; s < up_to; s++) {
        sum += chainPropagatorQStar[scalar_index(x,y,z,N - s)] * chainPropagatorQ[scalar_index(x,y,z,s)] * delta_s;
    }
    sum /= (compute_Q() * N);
    return sum;
}

double SCFT::compute_reduced_density_segment_B(int x, int y, int z) {
    double* chainPropagatorQ = q_propagator->get_chain_propagator();
    double* chainPropagatorQStar = q_dagger_propagator->get_chain_propagator();
    double sum = 0.0;
    const int up_to = lround(f * N);
    for(int s = up_to; s < N; s++) {
        sum += chainPropagatorQStar[scalar_index(x,y,z,N - s)] * chainPropagatorQ[scalar_index(x,y,z,s)] * delta_s;
    }
    sum /= (compute_Q() * N);
    return sum;
}


void SCFT::Determine_Density_Differences() {
    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                density_A[scalar_index(x,y,z)] = compute_reduced_density_segment_A(x,y,z);
                density_B[scalar_index(x,y,z)] = compute_reduced_density_segment_B(x,y,z);
            }
        }
    }
}


void SCFT::Run_Propagators() {
    s = 0;
    q_propagator->reset_time();
    q_dagger_propagator->reset_time();
    for(int i = 0; i < N_s; i++) {
        perform_timestep();
    }
}

void SCFT::Update_Fields() {
    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                double compressibility_condition = 0.5 * (w_A[scalar_index(x,y,z)] + w_B[scalar_index(x,y,z)] - chiN);
                double w_a_out = chiN * density_B[scalar_index(x,y,z)] + compressibility_condition;
                double w_b_out = chiN * density_A[scalar_index(x,y,z)] + compressibility_condition;
                double w_a = (w_A[scalar_index(x,y,z)] * (1.0 - mixing_parameter) + mixing_parameter * w_a_out);
                double w_b = (w_B[scalar_index(x,y,z)] * (1.0 - mixing_parameter) + mixing_parameter * w_b_out);
                w_A[scalar_index(x,y,z)] = w_a;
                w_B[scalar_index(x,y,z)] = w_b;
                //std::cout << "bob!" << '\n';
            }
        }
    }
    double w_a_b_average = 0.0;
    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                w_a_b_average += w_A[scalar_index(x,y,z)] + w_B[scalar_index(x,y,z)];
            }
        }
    }
    w_a_b_average /= (2.0 * NX * NY * NZ);
    for(int i = 0; i < NX * NY * NZ; i++) {
        w_A[i] -= w_a_b_average;
        w_B[i] -= w_a_b_average;
    }
}

double SCFT::Determine_Error() {
    double sum = 0.0;
    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                double compressibility_condiiton = 0.5 * (w_A[scalar_index(x,y,z)] + w_B[scalar_index(x,y,z)] - chiN);
                double w_a_out = chiN * density_B[scalar_index(x,y,z)] + compressibility_condiiton;
                double w_b_out = chiN * density_A[scalar_index(x,y,z)] + compressibility_condiiton;
                double difference_in_w_A = fabs(pow(w_a_out - w_A[scalar_index(x,y,z)], 2.0));
                double difference_in_w_B = fabs(pow(w_b_out - w_B[scalar_index(x,y,z)], 2.0));
                sum = sum + (difference_in_w_A + difference_in_w_B) * delta_x;
            }
        }
    }
    sum = sqrt(1.0 / ((NX * NY * NZ) * delta_x) * sum) / std::max(1.0, chiN);
    //output_fields.close();
    return sum;
}

void SCFT::Run() {
    if(field_type == "scft") {

        int index = 1.0;
        double field_error_threshold = pow(10.0,-3.0);
        double variance_threshold = pow(10.0,-5.0);
        double field_error = Determine_Error();
        double variance_error = sqrt(Determine_Variance_Total());
        while(field_error > field_error_threshold || variance_error > variance_threshold) {
            if(index == 1) {
                Determine_Density_Differences();
                output_lbm_data("output/0.csv", true, true);
            }
            double error_field = Determine_Error();
            double error_variance = Determine_Variance_Total();
            std::cout << index << " - Error: " << error_field << '\n';
            std::cout << index << " - Stdev: " << sqrt(error_variance) << '\n';
            Run_Propagators();
            Determine_Density_Differences();
            Update_Fields();
            std::string file_name = "output/" + std::to_string(index) + ".csv";
            if(index % 200 == 0 || index < 20) {
                output_lbm_data(file_name, true, true);
            }
            if(!stable()) {
                std::cout << "Warning: Numerical stability compromised |u_max|>=0.577" << '\n';
            }
            index++;
            if(index > 2000) {
                break;
            }
            field_error = Determine_Error();
            variance_error = sqrt(Determine_Variance_Total());
        }
        std::cout << index << " - Error: " << field_error << '\n';
        std::cout << index << " - Stdev: " << variance_error << '\n';
        output_lbm_data("output/output.csv", true, true);
    } else {
        Run_Propagators();
        output_lbm_data("output/output.csv", true, true);
    }

}

bool SCFT::stable() {
    return (q_propagator->stable() && q_dagger_propagator->stable());
}

bool SCFT::u_stable() {
    return (q_propagator->u_stable() && q_dagger_propagator->u_stable());
}

bool SCFT::relaxation_stable() {
    return q_propagator->relaxation_stable();
}

void SCFT::output_parameters() {
    q_propagator->output_parameters();
}

void SCFT::output_field_parameters() {
    q_propagator->output_field_parameters();
}

double SCFT::find_u_max() {
    return q_propagator->find_u_max();
}
