#include "csv_printer.hpp"

#include <utility>

CSVPrinter::CSVPrinter(std::string file_prefix, std::string file_extension) : 
    m_file_prefix{std::move(file_prefix)}, m_file_extension(std::move(file_extension))
{

}

void CSVPrinter::print(const std::vector<double> &x, const std::vector<double> &y, std::string file_name, unsigned num_ghost) const
{
    std::ofstream file;
    file.open(m_file_prefix + file_name + "." + m_file_extension);

    char delim = ',';

    // header
    file << "x" << delim << "y" << std::endl;

    if (x.size() != y.size()) {
        throw std::runtime_error("x and y need to be the same size");
        return;
    }

    for (size_t i = num_ghost; i < x.size() - num_ghost; ++i) {
        file << x[i] << delim << y[i] << std::endl; // trailing endl shoudn't cause problems
    }

    file.close();
}

void CSVPrinter::print_mat(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& m, std::string file_name, Eigen::Index num_ghost) const
{
    std::ofstream file;
    file.open(m_file_prefix + file_name + "." + m_file_extension);

    // we are only interested in the real part
    Eigen::MatrixXd mp(m.rows()-2*num_ghost, m.cols()-2*num_ghost);

    for (Eigen::Index j = num_ghost; j < m.cols() - num_ghost; ++j) {
        for (Eigen::Index i = num_ghost; i < m.rows() - num_ghost; ++i) {
            mp(i-num_ghost, j-num_ghost) = m(i, j).real();
        }
    }

    char delim = ',';

    for (Eigen::Index i = 0; i < mp.rows(); ++i) {
        for (Eigen::Index j = 0; j < mp.cols(); ++j) {
            file << mp(i, j);
            if (j != mp.cols() - 1) { // prevent trailing comma
                file << delim;
            }
        }
        file << std::endl;
    }
    
    file.close();
}