#include "genome.h"
#include <iostream>
#include <fstream>
#include <system_error>
#include <string>
#include <algorithm>


// class to read a DNA fasta file and concatenate its contents into a single string.



Genome::Genome(std::string arg, bool sequence) {
 // first case, arg represents actual sequence
    if (sequence) {
        std::transform(arg.begin(), arg.end(), arg.begin(), ::toupper); //make sure everything is upper case
        if (is_valid_dna(arg)) { // make sure all characters are 'G' 'C' 'A' 'T'
            genome = arg;
        } // now can proceed with other part of constructor.
    
    }
    else if (!sequence) { // for readability only
    
        std::string line;
        std::ifstream in_file {arg}; // note, according to RAII, this will close file if exception is thrown.
        if (!in_file.is_open()) { // if file is not found..
            throw std::system_error(errno, std::generic_category(), arg);  // this will terminate the program
        }
        else {
            while(std::getline(in_file, line)) {            
                if (line[0] != '>') {
                    genome += line;
                }
            
            }
        
        }
        if (in_file.bad()) { // any other problem reading file, terminate program.
            throw std::system_error(errno, std::generic_category(), arg);
        }
        in_file.close();
        std::transform(genome.begin(), genome.end(), genome.begin(), ::toupper);
        is_valid_dna(genome); // this will throw an exception if non DNA characters.
    
        make_complements();
    }
}
    
    
void Genome::make_complements(){
    complement = std::string(genome.size(), '\0');
    

    for (size_t i {0}; i < genome.size(); ++i){
        complement[i] = comp[genome[i]]; // make complement
    }
    rev_comp = complement; // preparation, this won't be the final result, see next line.
    std::reverse(rev_comp.begin(), rev_comp.end()); // reverse it
}


const std::string& Genome::get_genome() const{ // don't allow modification
    return genome;
    }
    
const std::string& Genome::get_rev_comp() const {
    return rev_comp;
}

const std::string& Genome::get_complement() const {
    return complement;
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
    
bool Genome::is_valid_dna(std::string& st) {
    if (!(std::all_of(st.begin(), st.end(), [](char x) {return (x == 'G' || x == 'C' || x == 'T' || x == 'A');}))) {
        throw std::runtime_error("file or string contains invalid characters");
        return false;
    }
    return true;
}
    
//}
    




