#include <stdio.h>
#include <iostream>
#include <fstream>
#include <streambuf>
//Now Linux only.
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "LBM.hpp"
#include "vector3.hpp"
//RapidJSON files.
#include "document.h"
#include "writer.h"
#include "stringbuffer.h"
#include "csv.h"
#include "../SCFT/SCFT.hpp"
#include <math.h>

using namespace rapidjson;

inline bool file_exists (const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

int main(int argc, char** argv) {
    if(!file_exists("options.json")) {
        std::cout << "Please ensure that options.json exists. If not, it can be obtained from root directory of GitHub repo." << '\n';
        return -1;
    }

    std::cout << "Do you want to clean the previous run? (1 - Yes, 0 - No): ";
    int choice;
    std::cin >> choice;
    if(choice == 1) {
        system("rm -rf output");
        system("mkdir output");
    }
    std::ifstream t("options.json");
    std::string str((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());
    rapidjson::Document d;
    d.Parse(str.c_str());
    Value& m_c_s = d["c_s"];
    Value& m_chain_length = d["chain_length"];
    Value& m_box_length_rg = d["box_length_rg"];
    Value& m_velocity_set = d["velocity_set"];
    Value& m_boundary_conditions = d["boundary_conditions"];
    Value& m_gamma_dot = d["gamma_dot"];
    Value& m_field_type = d["field_type"];
    Value& m_f = d["f"];
    Value& m_chiN = d["chiN"];
    Value& m_N_s = d["N_s"];
    auto grid_size = d["grid_size"].GetArray();
    int NX = grid_size[0].GetInt();
    int NY = grid_size[1].GetInt();
    int NZ = grid_size[2].GetInt();
    int chain_length = m_chain_length.GetInt();
    double c_s = m_c_s.GetDouble();
    double N_s = m_N_s.GetDouble();
    double box_length_rg = m_box_length_rg.GetDouble();
    const double R_g = NX / box_length_rg; // Using lattice units L = NX * delta_x, delta_x = 1, L = NX
    const double b = sqrt(6) / sqrt(chain_length) * R_g;
    //double tau_LB = 3 * b * b / 6.0 + 0.5;
    double tau_LB = b * b / (6.0 * c_s * c_s) + 0.5;
    double gamma_dot = m_gamma_dot.GetDouble();
    double viscosity = c_s * c_s * (tau_LB - 0.5);
    double radius_gyration = NX / box_length_rg;
    int scale = 1;
    int runs = chain_length * scale * scale * scale;
    double f = m_f.GetDouble();
    double chiN = m_chiN.GetDouble();
    std::string velocity_set = m_velocity_set.GetString();
    std::string boundary_conditions = m_boundary_conditions.GetString();
    std::string field_type = m_field_type.GetString();
    if(velocity_set != "D3Q15" && velocity_set != "D3Q27" && velocity_set != "D2Q9" && velocity_set != "D3Q7") {
        std::cout << "Error: Please specify a valid velocity set such as D3Q15,D3Q27,D3Q7 or D2Q9." << '\n';
        return -1;
    }
    if(boundary_conditions != "lees_edwards" && boundary_conditions != "periodic" && boundary_conditions != "couette") {
        std::cout << "Errors: boundary_conditions in options.json can either be: periodic, Couette (D2Q9 only) or lees_edwards (Lees-Edwards Shear, Please see research paper by Alexander Wagner)";
        return -1;
    }
    std::cout << "Shear rate (gamma_dot): " << gamma_dot << '\n';
    if(NZ != 1 && velocity_set == "D2Q9") {
        std::cout << "Warning: NZ=1 for D2Q9.";
        return -1;
    }
    SCFT *solver = new SCFT(
            NX,NY,NZ,
            velocity_set,c_s,boundary_conditions,
            gamma_dot,chain_length,N_s,f,chiN,box_length_rg,
            field_type,0.1
            );
    /*for(int i = 0; i < argc; i++) {
        if(std::string(argv[i]) == "generate_ic") {
            solver->output_lbm_data("ic.csv", true);
            std::cout << "Generated ic.csv" << '\n';
            return 0;
        }
    }
    if(file_exists("ic.csv")) {
        std::cout << "Loading initial conditions" << '\n';
        io::CSVReader<4> in("ic.csv");
        in.read_header(io::ignore_extra_column, "p","u_x","u_y","u_z");
        double density,u_x,u_y,u_z;
        for(int i = 0; i < NX; i++) {
            for(int j = 0; j < NY; j++) {
                for(int k = 0; k < NZ; k++) {
                    in.read_row(density,u_x,u_y,u_z);
                    solver->set_density(i,j,k,density);
                    solver->set_velocity(i,j,k,u_x,u_y,u_z);
                }
            }
        }
        std::cout << "Loaded initial conditions" << '\n';
    } else {
        std::cout << "Using default of p=1 for all x,y,z and u(x,t=0)=0 for all x,y,z. (Steady state)" << '\n';
        std::cout << "If you wish to use your own initial conditions, please run the program but with command: generate_ic as a argument which will output ic.csv in format of p,u_x,u_y,u_z, assume indexes are incrementing i,j,k for i<NX,j<NY and k<NZ" << '\n';
    }*/
    //Equation 3.5 with delta t = 1, LBM Principles and Practice book.
    std::cout << "----Lattice Boltzmann parameters----" << '\n';
    solver->output_parameters();
    //Equation 4.4.49
    std::cout << "----------Field parameters----------" << '\n';
    std::cout << "Number of integration steps N_s: " << N_s << '\n';
    solver->output_field_parameters();
    solver->output_lbm_data("output/0.csv", true, true);
    solver->Run();
    solver->output_max();
    solver->output_min();
    std::cout << "-----After simulation------" << '\n';
    if(solver->stable()) {
        std::cout << "Stable |u_max|<0.577 for D2Q9,D3Q15 and D3Q27 and relaxation>0.5: " << "Yes" << '\n';
    } else if(!solver->u_stable()) {
        std::cout << "Stable |u_max|<0.577 for D2Q9,D3Q15 and D3Q27 and relaxation>0.5: " << "No (|u_max|=" << solver->find_u_max() << ")" << '\n';
    } else if(!solver->relaxation_stable()) {
        std::cout << "Stable |u_max|<0.577 for D2Q9,D3Q15 and D3Q27 and relaxation>0.5: " << "No (relaxation=" << solver->tau_LB << ")" << '\n';
    }
    std::cout << std::endl;
    return 0;
}
