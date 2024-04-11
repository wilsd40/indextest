#ifndef GENOME_H
#define GENOME_H

#include <iostream>
#include <fstream>
#include <system_error>
#include <string>
#include <map>


class Genome{
public:
    std::map<char, char> comp = {{'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'}} ;
    Genome(std::string arg, bool sequence = true);
    const std::string& get_genome() const;
    const std::string& get_rev_comp() const;
    const std::string& get_complement() const;
    void make_complements();
    char get_pos(int);
    std::string_view get_genome(int, size_t);

    int get_size();
    bool is_valid_dna(std::string&);
    std::string genome;
    std::string complement;
    std::string rev_comp;




};

#endif // GENOME_H
