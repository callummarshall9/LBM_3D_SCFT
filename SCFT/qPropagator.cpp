//
// Created by callummarshall on 25/11/2019.
//

#include "qPropagator.hpp"

qPropagator::qPropagator(int nx, int ny, int nz, std::string velocity_set, double c_s, std::string boundary_condition,
        double gamma_dot, double N, int N_s, double box_length_rg, std::string field_type, double f, double chiN,
        double* w_A, double* w_B):
    LBM(nx,ny,nz,velocity_set,c_s,boundary_condition,gamma_dot,N,N_s,box_length_rg,field_type,f,chiN,w_A,w_B) {
}

double qPropagator::field(int x, int y, int z) {
    //Reference [1] - Fields omega_{pm} = N lowercase omega_{pm}
    if(field_type == "scft") {
        if(s < N * f) {
            return w_A[scalar_index(x,y,z)] / N;
        } else {
            return w_B[scalar_index(x,y,z)] / N;
        }
    } else {
        return LBM::field(x,y,z);
    }

}
