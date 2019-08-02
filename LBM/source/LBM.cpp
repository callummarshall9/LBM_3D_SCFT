#include "LBM.hpp"
#include "vector3.cpp"
#include <iomanip>
#include <cmath>
#include <assert.h>

LBM::LBM(int grid_size, std::string velocity_set, double m_c_s, double m_tau, std::string m_boundary_condition, double m_gamma_dot) : NX(grid_size),NY(grid_size),NZ(grid_size),
        c_s(m_c_s), tau_LB(m_tau), boundary_condition(m_boundary_condition), gamma_dot(m_gamma_dot), velocity_set(velocity_set) {
    set_velocity_set(velocity_set);
    initialise();
}

LBM::LBM(int nx, int ny, int nz, std::string velocity_set, double m_c_s, double m_tau, std::string m_boundary_condition, double m_gamma_dot) : NX(nx), NY(ny), NZ(nz),
        c_s(m_c_s), tau_LB(m_tau), boundary_condition(m_boundary_condition), gamma_dot(m_gamma_dot), velocity_set(velocity_set) {
    set_velocity_set(velocity_set);
    initialise();
}

void LBM::initialise() {
    int box_flatten_length = NX * NY * NZ;
    int distributions_flatten_length = box_flatten_length * direction_size;
    density_field = new double[box_flatten_length];
    velocity_field = new vector3<double>[box_flatten_length];
    particle_distributions_star = new double[distributions_flatten_length];
    particle_distributions = new double[distributions_flatten_length];
    field_plus = new double[box_flatten_length];
    field_minus = new double[box_flatten_length];
    for(int i = 0; i < NX * NY * NZ; i++) {
        density_field[i] = 1.0;
    }
    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                for(int i = 0; i < direction_size; i++) {
                    particle_distributions_star[scalar_index(x,y,z,i)] = weights[i] * density_field[scalar_index(x,y,z)];
                    particle_distributions[scalar_index(x,y,z,i)] = weights[i] * density_field[scalar_index(x,y,z)];
                }
            }
        }
    }
}

double LBM::calculate_feq(int x, int y, int z, int i) {
    double dot_product = (double)velocity_field[scalar_index(x,y,z)].x * (double)directions[i].x + (double)velocity_field[scalar_index(x,y,z)].y * (double)directions[i].y +
                         (double)velocity_field[scalar_index(x,y,z)].z * (double)directions[i].z;
    //Equation 3.4 with c_s^2 = 1/3
    //double feq = weights[i] * density_field[scalar_index(x,y,z)] * (1.0 + dot_product / (c_s * c_s) + dot_product * dot_product / (2 * c_s * c_s * c_s * c_s) - velocity_field[scalar_index(x,y,z)].norm_square() / (2 * c_s * c_s));
    //Equation (21) from Lattice Boltzmann Method for multiscale self-consistent field theory
    //simulations of block copolymers by Hsieh Chen, YongJoo Kim and Alfredo Alexander-Katz.
    double feq = weights[i] * density_field[scalar_index(x,y,z)];
    return feq;
}

double LBM::calculate_feq(int x, int y, int z, int i, double u_le_x) {
    double dot_product = (velocity_field[scalar_index(x,y,z)].x + u_le_x) * directions[i].x + velocity_field[scalar_index(x,y,z)].y * directions[i].y +
                         velocity_field[scalar_index(x,y,z)].z * directions[i].z;
    double norm_square = (velocity_field[scalar_index(x,y,z)].x + u_le_x) * (velocity_field[scalar_index(x,y,z)].x + u_le_x) + velocity_field[scalar_index(x,y,z)].y * directions[i].y +
                            velocity_field[scalar_index(x,y,z)].z * directions[i].z;
    //Equation 3.4 with c_s^2 = 1/3
    double feq = weights[i] * density_field[scalar_index(x,y,z)] * (1.0 + dot_product / (c_s * c_s) + dot_product * dot_product / (2 * c_s * c_s * c_s * c_s) - norm_square / (2 * c_s * c_s));
    return feq;
}

void LBM::set_velocity_set(std::string velocity_set) {
    if(velocity_set == "D3Q15") {
        direction_size = 15;
        directions = new vector3<int>[15]{ {0,0,0},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,1},{-1,-1,-1},{1,1,-1},{-1,-1,1},{1,-1,1},{-1,1,-1},{-1,1,1},{1,-1,-1} };
        weights = new double[15] { 2.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0  };
    } else if(velocity_set == "D3Q27") {
        direction_size = 27;
        directions = new vector3<int>[27]{{0,0,0}, {1,0,0},{-1,0,0},
        {0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,-1,0},{1,0,1},
        {-1,0,-1},{0,1,1},{0,-1,-1},{1,-1,0},{-1,1,0},{1,0,-1},{-1,0,1},
        {0,1,-1},{0,-1,1},{1,1,1},{-1,-1,-1},{1,1,-1},{-1,-1,1},{1,-1,1},
        {-1,1,-1},{-1,1,1},{1,-1,-1}};
        weights = new double[27] {8.0/27.0,2.0/27.0,2.0/27.0,2.0/27.0,
            2.0/27.0,2.0/27.0,2.0/27.0, 1.0/54.0,1.0/54.0,1.0/54.0
        ,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,
        1.0/54.0,1.0/54.0,1.0/54.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0
        ,1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0};
    } else if(velocity_set == "D2Q9") {
        direction_size = 9;
        directions = new vector3<int>[9]{{1,0,0},{0,1,0},{-1,0,0},{0,-1,0},{1,1,0},{-1,1,0},{-1,-1,0},{1,-1,0},{0,0,0}};
        weights = new double[9]{1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,4.0/9.0};
    }
    lookup_reverse();

}

void LBM::lookup_reverse() {
    reverse_indexes = new int[direction_size];
    for(int i = 0; i < direction_size; i++) {
        for(int j = 0; j < direction_size; j++) {
            if(directions[i].x == -directions[j].x && directions[i].y == -directions[j].y && directions[i].z == -directions[j].z) {
                reverse_indexes[i] = j;
            }
        }
    }
}

LBM::~LBM() {
    delete[] density_field;
    delete[] particle_distributions;
}

void LBM::output_array(double* array) {
    std::cout << "x,y,z value" << std::endl;
    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                std::cout << x << "," << y << "," << z << ": " << array[this->scalar_index(x,y,z)] << std::endl;
            }
        }
    }
}

void LBM::output_density() {
    this->output_array(this->density_field);
}

void LBM::output_velocity() {
    int z_index = 0;
    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            std::cout << velocity_field[scalar_index(x,y,z_index)].x << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void LBM::set_velocity(int x_field, int y_field, int z_field, double u_x, double u_y, double u_z) {
    velocity_field[scalar_index(x_field,y_field,z_field)].x = u_x;
    velocity_field[scalar_index(x_field,y_field,z_field)].y = u_y;
    velocity_field[scalar_index(x_field,y_field,z_field)].z = u_z;
}

void LBM::set_density(int x_field, int y_field, int z_field, double density) {
    density_field[scalar_index(x_field,y_field,z_field)] = density;
}


void LBM::compute_density_momentum_moment() {
    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                //Equation 3.1
                //Equation (19) from Lattice Boltzmann Method for multiscale self-consistent field theory
                //simulations of block copolymers by Hsieh Chen, YongJoo Kim and Alfredo Alexander-Katz.
                double new_density = 0;
                vector3<double> u;
                for(int i = 0; i < direction_size; i++) {
                    new_density += particle_distributions[scalar_index(x,y,z,i)];
                    u.x += particle_distributions[scalar_index(x,y,z,i)] * directions[i].x;
                    u.y += particle_distributions[scalar_index(x,y,z,i)] * directions[i].y;
                    u.z += particle_distributions[scalar_index(x,y,z,i)] * directions[i].z;
                }
                density_field[scalar_index(x,y,z)] = new_density;
                velocity_field[scalar_index(x,y,z)].x = u.x / new_density;
                velocity_field[scalar_index(x,y,z)].y = u.y / new_density;
                velocity_field[scalar_index(x,y,z)].z = u.z / new_density;
            }
        }
    }
}

void LBM::stream() {
    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                for(int i = 0; i < direction_size; i++) {
                    //Periodic boundary conditions taken from Taylor green in Chapter 13.
                    if (boundary_condition == "periodic") {
                        //X position Minus the Direction (xmd) applies to y and z.
                        int xmd = (NX + x - (int)directions[i].x) % NX;
                        int ymd = (NY + y - (int)directions[i].y) % NY;
                        int zmd = (NZ + z - (int)directions[i].z) % NZ;
                        particle_distributions[scalar_index(x,y,z,i)] = particle_distributions_star[scalar_index(xmd,ymd,zmd,i)];
                        //Equation 3.10 with periodic boundary conditions.
                    } else if(boundary_condition == "couette") {

                        //Using Bounce back approach for top & bottom walls..
                        if (y==0 && directions[i].y == 1) {
                            //Bottom Wall.
                            //Equation 5.27 from LBM Principles and Practice.
                            particle_distributions[scalar_index(x,y,z,i)]=particle_distributions_star[scalar_index(x,y,z,reverse_indexes[i])];
                        } else if (y==NY-1 && directions[i].y == -1) {
                            //Top wall
                            //Equation 5.28 from LBM Principles and Practice.
                            //coefficients of Equation 5.28 calculated from footnote 17.
                            double u_max = 0.1;
                            particle_distributions[scalar_index(x,y,z,i)]=particle_distributions_star[scalar_index(x,y,z,reverse_indexes[i])]+directions[i].x * 2 * weights[i] / (c_s * c_s) * u_max;
                            //particle_distributions[scalar_index(x,y,z,(4-1))]=previous_particle_distributions[scalar_index(x,y,z,(2-1))];
                            //particle_distributions[scalar_index(x,y,z,(7-1))]=previous_particle_distributions[scalar_index(x,y,z,(5-1))]-(1.0/6.0)*u_max;
                            //particle_distributions[scalar_index(x,y,z,(8-1))]=previous_particle_distributions[scalar_index(x,y,z,(6-1))]+(1.0/6.0)*u_max;
                        } else {
                            //Chapter 13 Taylor-Green periodicity from same book.
                            int xmd = (NX + x - directions[i].x) % NX;

                            //int ymd = (NY + y - (int)directions[i].y) % NY;
                            int ymd = y - directions[i].y;
                            int zmd = (NZ + z - directions[i].z) % NZ;
                            particle_distributions[scalar_index(x,y,z,i)] = particle_distributions_star[scalar_index(xmd,ymd,zmd,i)];
                        }
                    } else if(boundary_condition == "lees_edwards") {
                            double d_x = gamma_dot * (double)NY * (double)time;
                            int d_x_I = d_x;
                            double d_x_R = d_x - d_x_I;
                            d_x_I = d_x_I % NX;
                            double s_1 = d_x_R;
                            double s_2 = 1.0 - d_x_R;
                            //s_1=0;s_2=1;
                            int ymd = (NY + y - directions[i].y) % NY;
                            int zmd = (NZ + z - directions[i].z) % NZ;
                            if(y==0 && directions[i].y == 1) {
                                //Bottom Wall.
                                //Equation (17) from Less-Edwards boundary conditions for lattice Boltzmann suspension simulations
                                //by Eric Lorenz and Alfons G. Hoekstra
                                    int x_pos;
                                    x_pos = (x+d_x_I + 2 * NX - directions[i].x) % NX;
                                    int x_shifted;
                                    x_shifted = (x+d_x_I + 2 * NX +1 - directions[i].x) % NX;
                                    double galilean_transformation_pos = particle_distributions_star[scalar_index(x_pos,ymd,zmd,i)] + calculate_feq(x_pos,ymd,zmd,i,-1 * gamma_dot * (double)NY) - calculate_feq(x_pos,ymd,zmd,i,0);
                                    double galilean_transformation_shift = particle_distributions_star[scalar_index(x_shifted,ymd,zmd,i)] + calculate_feq(x_shifted,ymd,zmd,i,-1 * gamma_dot * (double)NY) - calculate_feq(x_shifted,ymd,zmd,i,0);
                                    //Equation (18) from same paper.
                                    particle_distributions[scalar_index(x,y,z,i)]= s_1 * galilean_transformation_shift +
                                                                                   s_2 * galilean_transformation_pos;
                                } else if(y==NY-1 && directions[i].y == -1) {
                                    //Equation (17) from Less-Edwards boundary conditions for lattice Boltzmann suspension simulations
                                    //by Eric Lorenz and Alfons G. Hoekstra
                                    //Top Wall.
                                    int x_pos;
                                    x_pos = (x-d_x_I + 2 * NX - directions[i].x) % NX;
                                    int x_shifted;
                                    x_shifted = (x-d_x_I + 2 * NX -1 - directions[i].x) % NX;
                                    double galilean_transformation_pos = particle_distributions_star[scalar_index(x_pos,ymd,zmd,i)] + calculate_feq(x_pos,ymd,zmd,i,1 * gamma_dot * (double)NY) - calculate_feq(x_pos,ymd,zmd,i,0);
                                    double galilean_transformation_shift = particle_distributions_star[scalar_index(x_shifted,ymd,zmd,i)] + calculate_feq(x_shifted,ymd,zmd,i,1 * gamma_dot * (double)NY) - calculate_feq(x_shifted,ymd,zmd,i,0);
                                    //Equation (18) from same paper.
                                    particle_distributions[scalar_index(x,y,z,i)]= s_1 * galilean_transformation_shift +
                                                                                   s_2 * galilean_transformation_pos;
                                } else {
                                    //Interior points. (Enforce periodicity along x axis);.
                                    int xmd = (NX + x - directions[i].x) % NX;
                                    particle_distributions[scalar_index(x,y,z,i)] = particle_distributions_star[scalar_index(xmd,ymd,zmd,i)];
                                }


                    }
                }
            }
        }
    }
}

void LBM::output_indices_file() {
    std::ofstream output("output/indices.csv");
    output << "x,y,z" << '\n';
    for(int i = 0; i < NX; i++) {
        for(int j = 0; j < NY; j++) {
            for(int k = 0; k < NZ; k++) {
                output << i << "," << j << "," << k << '\n';
            }
        }
    }
    std::cout << std::endl;
    output.close();
}

double sech(double x) {
    return 1.0 / cosh(x);
}

void LBM::set_field_paramaeters(double N, double Length_Rg, std::string m_field_type) {
    chain_length = N;
    box_length_rg = Length_Rg;
    const double delta_x = 1.0;//Lattice units.
    const double L = (double)NX * delta_x;
    R_g = L / box_length_rg;
    const double b = sqrt(6) / sqrt(chain_length) * R_g;
    tau_LB = (3.0 * b * b / 6.0 + 0.5);
    field_type = m_field_type;
}

double LBM::field(int x, int y, int z) {
    const double delta_x = 1.0;
    const double L = (double)NX * delta_x;
    double x_value = x * delta_x;
    double midpoint = NX / (2.0);
    double field;
    if(field_type == "fixed_sine") {
        double A = 10.0;
        double M = 1.0;
        field = A * sin(M * 2 * M_PI * (x_value - midpoint) / L) / chain_length;
    } else if(field_type == "fixed_sech") {
        field = (1.0  - 2.0 * pow(sech(3.0 * (x_value - midpoint) / (2 * R_g)),2.0)) / chain_length;
    }
    return field;
}

void LBM::collision() {//Performs the collision step.
    const double tauinv = 1.0 / tau_LB;
    const double omtauinv = 1.0-tauinv;

    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                for(int i = 0; i < direction_size; i++) {

                    double R = -field(x,y,z) * density_field[scalar_index(x,y,z)];
                    double feq = calculate_feq(x,y,z,i);
                    //Equation 20
                    particle_distributions_star[scalar_index(x,y,z,i)] = particle_distributions[scalar_index(x,y,z,i)] +  tauinv * (feq - particle_distributions[scalar_index(x,y,z,i)])  + weights[i] * R;
                    //Equation 3.9
                    //previous_particle_distributions[scalar_index(x,y,z,i)] = omtauinv * particle_distributions[scalar_index(x,y,z,i)] + tauinv * feq + weights[i] * R;
                }
            }
        }
    }
}

double gauss_noise(double x) {
    double mean = 0;
    double stdev = 1.0;
    return 1.0 / (stdev * sqrt(2 * M_PI)) * exp(-1.0 * pow(x - mean,2.0) / (2 * stdev * stdev));
}

void LBM::perform_field_calculations() {
    double C = 1.0;
    double T = 1.0;
    double florry_higgs = 0.5;
    double N = 1000;//Total number of statistical segments.
    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                double density_segment_A = 0;
                double density_segment_B = 0;
                double density_minus = density_segment_A - density_segment_B;
                field_minus[scalar_index(x,y,z)] = field_minus[scalar_index(x,y,z)] - T * C * (-1.0 * density_minus + 2.0 / (florry_higgs * N) * field_minus[scalar_index(x,y,z)]) + gauss_noise(time);
            }
        }
    }
}

int LBM::get_time() {return time;}

void LBM::perform_timestep() {
    time++;
    compute_density_momentum_moment();
    collision();//particle_distributions=previous_
    stream();//previous=particle_
}

void LBM::output_f_array(double* f_array, int z_index) {
    for(int i = 0; i < direction_size; i++) {
        std::cout << "ans(:,:," << z_index << "," << (i+1) << ")" << '\n';
        std::cout << '\n';
        for(int x = 0; x < NX; x++) {
            for(int y = 0; y < NY; y++) {
                std::cout << f_array[scalar_index(x,y,z_index,i)] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << '\n';
    }
    std::cout << std::endl;
}


void LBM::output_lbm_data(std::string filename, bool header, bool output_indices) {
    std::ofstream output_stream;
  output_stream.open (filename, std::ofstream::out | std::ofstream::app);
    if(header) {
        if(output_indices) {
            output_stream << "x,y,z,";
        }
        output_stream << "p,u_x,u_y,u_z" << '\n';
    }
    double a;
    double b;
    for(int x = 0; x < NX; x++) {
        for(int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++) {
                if(output_indices) {
                    output_stream << x << "," << y << "," << z << ",";
                }
                output_stream << density_field[this->scalar_index(x,y,z)] << "," <<
                    velocity_field[this->scalar_index(x,y,z)].x << "," <<
                    velocity_field[this->scalar_index(x,y,z)].y << "," <<
                    velocity_field[this->scalar_index(x,y,z)].z << "," << field(x, a, b) << '\n';
            }
        }
    }
    output_stream.close();
}
