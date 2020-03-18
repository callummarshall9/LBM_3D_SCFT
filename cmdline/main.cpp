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
#include "../field_types.cpp"
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
    double mixing_parameter = d["mixing_parameter"].GetDouble();
    double field_error_threshold = d["field_error_threshold"].GetDouble();
    double variance_threshold = d["variance_threshold"].GetDouble();
    auto grid_size = d["grid_size"].GetArray();
    int NX = grid_size[0].GetInt();
    int NY = grid_size[1].GetInt();
    int NZ = grid_size[2].GetInt();
    double chain_length = d["chain_length"].GetDouble();
    double c_s = d["c_s"].GetDouble();
    int N_s = d["N_s"].GetInt();
    double box_length_rg = d["box_length_rg"].GetDouble();
    double gamma_dot = d["gamma_dot"].GetDouble();
    int scale = 1;
    double f = d["f"].GetDouble();
    double chiN = d["chiN"].GetDouble();
    std::string velocity_set = d["velocity_set"].GetString();
    std::string boundary_conditions = d["boundary_conditions"].GetString();
    std::string m_field_type = d["field_type"].GetString();
    std::string field_initial = d["field_initial"].GetString();
    field_types field_type;
    if(m_field_type == "fixed_sine") {
        field_type = field_types::fixed_sine;
    } else if(m_field_type == "fixed_sech") {
        field_type = field_types::fixed_sech;
    } else if(m_field_type == "scft") {
        field_type = field_types::scft;
    } else {
        std::cout << "Error: Please use field_type value of fixed_sine,field_sech or scft." << '\n';
        return -1;
    }
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
            field_type,mixing_parameter,field_error_threshold, variance_threshold, field_initial
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
