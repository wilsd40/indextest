#ifndef GENOME_H
#define GENOME_H

#include <iostream>
#include <fstream>
#include <system_error>
#include <string>


class Genome{
public:
    Genome(std::string);
    const std::string& get_genome() const;
    char get_pos(int);
    std::string_view get_genome(int, size_t);
    int get_size();
private:
    std::string genome;
};

#endif // GENOME_H
