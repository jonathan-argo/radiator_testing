#include "gen.h"


void gen::checkSVD(const Eigen::MatrixXd& A) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A);
    int rank = svd.rank();
    std::cout << "Rank (SVD): " << rank << std::endl;
    if (rank == A.cols()) {
        std::cout << "Matrix is full rank (independent)." << std::endl;
    } else {
        std::cout << "Matrix is not full rank (dependent)." << std::endl;
    }
}