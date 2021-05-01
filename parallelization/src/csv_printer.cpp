#include "csv_printer.hpp"

#include <utility>

CSVPrinter::CSVPrinter(std::string file_prefix, std::string file_extension) : 
    m_file_prefix{std::move(file_prefix)}, m_file_extension(std::move(file_extension))
{

}

void CSVPrinter::print(const list& x, const list &y, std::string file_name, unsigned num_ghost) const
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

    for (ind i = num_ghost; i < x.size() - num_ghost; ++i) {
        file << x[i] << delim << y[i] << std::endl; // trailing endl shoudn't cause problems
    }

    file.close();
}