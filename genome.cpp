#include "genome.h"
#include <iostream>
#include <fstream>
#include <system_error>
#include <string>


// class to read a DNA fasta file and concatenate its contents into a single string.



Genome::Genome(std::string filename) { // in constructor, specify filename
    std::string line;
    std::ifstream in_file {filename}; // note, according to RAII, this will close file if exception is thrown.
    if (!in_file.is_open()) { // if file is not found..
        throw std::system_error(errno, std::generic_category(), filename);  // this will terminate the program
    }
    else {
        while(std::getline(in_file, line)) {            
            if (line[0] != '>') {
                genome += line;
            }
            
        }
        
    }
    if (in_file.bad()) { // any other problem reading file, terminate program.
        throw std::system_error(errno, std::generic_category(), filename);
    }
    in_file.close();
    
    }
const std::string& Genome::get_genome() const{ // don't allow modification
    return genome;
    }
    
char Genome::get_pos(int pos) {
    return genome.at(pos);
}
    
std::string_view Genome::get_genome(int start, size_t end){
    if (start < 0 || end >= genome.size()){
        return "";
    }
    return std::string_view(genome).substr(start, (end - start + 1));
    }


int Genome::get_size() {
    return genome.size();
    }
    




