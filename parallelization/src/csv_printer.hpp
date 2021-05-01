#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "commons.hpp"

class CSVPrinter 
{
private:
    std::string m_file_prefix;
    std::string m_file_extension;

    // static void printVector();
    // void print_list(std::ofstream& stream, const list& array, unsigned num_ghost);

public:
    CSVPrinter() : CSVPrinter("out/nsg_", "dat") {};
    CSVPrinter(std::string file_prefix, std::string extension);
    // void print(const State &state, double time) override;

    void print(const list& x, const list &y, std::string file_name, unsigned num_ghost=0) const;
};
