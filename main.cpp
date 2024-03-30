#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>
#include "genome.h"
#include <climits>
#include <chrono>
#include <numeric>
#include <optional>
#include <utility>
#include <cassert>



class Timer {
private:
    using Clock = std::chrono::steady_clock;
    using Second = std::chrono::duration<double, std::ratio<1>>;
    std::chrono::time_point<Clock> m_beg { Clock::now() };

public:
    void reset() {
        m_beg = Clock::now();
    }

    double elapsed() const {
        return std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - m_beg).count();
    }
};



class GenomeIndex{
    
    // note: ksize of 5 is a good compromise between speed of index creation (almost instantaneous) and search.
    // at least for a bacterial-sized search string.
    
using string_map = std::unordered_map<std::string, std::vector<int>>;
// this map provides constant access time, kmers are the string, and positions are in the vector.
using result_vector = std::vector<std::pair<int, int>>;
// results for query: [start position, length of identity].


private:
    string_map mp;

    
public:
    int ksize;
    
GenomeIndex(const std::string& te1, int ksize1) {
    make_index(te1, ksize1);
    this->ksize = ksize1;
    }
    
GenomeIndex(const Genome& ge, int ksize1){
    make_index(ge.get_genome(), ksize1);
    this->ksize = ksize1;
}
    
void make_index(const std::string& te, int ksize) {
    for (size_t i {0}; i < te.size() - ksize + 1; ++i) {
        std::string kmer = te.substr(i, ksize);
        mp[kmer].push_back(i); // no need the create an entry if one doesn't exist, it is already initialized
    }
    
}


const string_map& get_index() const{
    return mp;
}

void display() {
    if (mp.size() < 1000) {
    for (const auto& pair:mp){
        std::cout << pair.first << " :  [ ";
        for (const auto& pos:pair.second) {
            std::cout << pos << " , ";
        }
        std::cout << " ]" << std::endl;
    }
    }
    else {
        std::cout << "Structure too large to print" << std::endl;
    }
    
    
}
const int get_ksize() const{
    return this->ksize;
}
};

class GenomeIndexSearch {
using string_map = std::unordered_map<std::string, std::vector<int>>;
// this map provides constant access time, kmers are the string, and positions are in the vector.
using result_vector = std::vector<std::pair<int, int>>;
// results for query: [start position, length of identity].    

private:
    const std::string& ge;
    const string_map& gi;
    std::string quer;
    int ksize;
    int min_match;
    int direction; // 0, fwd, 1, rev, 2, don't know
    result_vector fwd_result;
    result_vector rev_result;
    result_vector result_return;


    
public:
   
    std::unordered_map<char, char> comp;
     
    GenomeIndexSearch(Genome& gearg, const GenomeIndex& giarg, std::optional<int> arg3 = std::nullopt): ge {gearg.get_genome()}, gi {giarg.get_index()},
        ksize {giarg.get_ksize()}{
        std::cout << "ksize = "<<ksize << std::endl;
        // constructor for using Genome object
        comp = {{'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'}} ; // when refactor, use this.
        
        
        // for reverse complement

        if (arg3.has_value()) {
            if (arg3 > 0 && arg3 < 3) {
                direction = arg3.value();
            }
        }
    }
    GenomeIndexSearch(std::string& gearg, const GenomeIndex& giarg, std::optional<int> arg3 = std::nullopt): ge {gearg}, gi {giarg.get_index()},
        ksize {giarg.get_ksize()}{
            std::cout << "ksize = "<<ksize << std::endl;
        // constructor for using Genome object
        comp = {{'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'}} ;
        
        
        // for reverse complement

        if (arg3.has_value()) {
            if (arg3 > 0 && arg3 < 3) {
                direction = arg3.value();
            }
        }
    }
// workhorse function for finding seed strings in index.  returns a vector of positions for which the string can be found.    
    std::vector<int> query_index(std::string query) {
        // make sure the check that query was length of index.
        std::string found;
        std::vector<int>rv;
        auto it1 = gi.find(query);
        if (it1 == gi.end()) {
            rv.push_back(-1);
        }
        else {
            found = it1->first;
            //std::cout << found << std::endl;
            const auto& ve {it1->second}; // ve will be vector from value part of key/value
            
            for (auto i:ve) { //iterate through hits
                //std::cout << i << std::endl;
                rv.push_back(i);
            }
        }
    return rv;
    }

// function to find primer in template.  calls query_index
bool search_string(std::string query, int dir, int min_match) {
    // this method will be called by another method to specify dir.  
    // search the subject with the query 
    // dir : 0 for forward primer (will start at 3' end and work backward on fwd strand
    //     : 1 for reverse primer (will make reverse complement, start at 5' end and work forward on fwd strand
    int end_nuc {-1};
    
    //need to check min_match vs query.
    std::cout << (int)query.length() << ": ksize :" <<ksize << std::endl;
    assert ((int)query.length() >= ksize);
    assert (min_match >= ksize);
    if (dir == 0) { // forward primer
        fwd_result.clear();
        //std::cout << (int)query.length() << std::endl;
        int start_nuc {(int)query.length() - ksize}; // 
        
        
        //std::cout << "Start nuc: "<<start_nuc<<std::endl;
 
        

        std::string qu_nucs = query.substr(start_nuc, ksize);
        //std::cout << "Will start search with: "<< qu_nucs << std::endl;

        auto indexstart = query_index(qu_nucs);
        if (indexstart.size() == 1 && indexstart.at(0) == -1) {
            std::cout << "search string not found" << std::endl;
            fwd_result.push_back(std::pair<int, int>{-1, -1});
            return 1;
            
        }
        for (auto fi:indexstart) {
        //auto fi = starttry.at(0);
        
        //for (auto fi:query_res){
            
            
            std::cout << "Starting at position "<<fi << std::endl;
           
        
   
            std::string result_string;
            int xtend {0};
            //std::cout << "conditions" << (start_nuc - xtend - 1) << " " << (fi - xtend - 1) << std::endl;
            // shouldn't have to check 3' bounds
            //bool found = false;
            while (((start_nuc - xtend - 1) >= 0 ) && (fi - xtend - 1) >=0){
                //std::cout << "------" << std::endl;
                std::cout << "before: bounds checks : start_nuc-xtend-1: " << start_nuc-xtend-1 << " fi - xtend - 1:  " << (fi - xtend - 1)<<std::endl;
                std::cout << xtend << " comparing " << query.at(start_nuc - xtend - 1) << " to " << ge.at(fi - xtend - 1) << std::endl;
                if ((query.at(start_nuc - xtend - 1)) == ge.at(fi - xtend - 1)) {
                    ++xtend;
                    //std::cout << "after: bounds checks : start_nuc-xtend-1: " << start_nuc-xtend-1 << " fi - xtend - 1:  " << (fi - xtend - 1)<<std::endl;
                }
                else {
                    break;
                }
            }

            std::string found = ge.substr(fi - xtend, xtend + ksize);
        
            //std::cout << "=====" << std::endl;
            if (xtend + ksize>= (min_match)){
                //std::cout << "Found: "<<found << std::endl;
                fwd_result.push_back(std::pair<int, int>(fi - xtend, xtend)); // need to think about this
            }
            }
    }
    
return 0;

}

bool check_fwd() {
    if (fwd_result.size() == 1 && fwd_result.at(0) == std::pair<int, int>(-1, -1)){
        return true;
    } else {
        return false;
    }
    
    
}

result_vector get_fwd() {
    return fwd_result;
}
    
    
void display_fwd_hits() {
    if (!check_fwd()){
        for (auto i: fwd_result) {
            int st = i.first;
            int xt = i.second;
            std::cout << "Position: " << st << " identity "<<xt<<"  Found "<<ge.substr(st, ksize+xt) << std::endl;
        }
    }
    }


};

std::vector<size_t> alt_find(std::string str, std::string sub) {
    std::vector<size_t> positions; // holds all the positions that sub occurs within str
    size_t pos = str.find(sub, 0);
    while(pos != std::string::npos)
    {
        positions.push_back(pos);
        pos = str.find(sub,pos+1);
    }
    return positions;
}

int compare_hits(std::vector<size_t> alt, GenomeIndexSearch obj){
    auto ve {obj.get_fwd()};
    int tot_mismatch {0};
    if (alt.size() == ve.size()){
        for (size_t i {0}; i < alt.size(); ++i) {
            std::cout << i << " " << ((int)alt.at(i)) << " : " << (ve.at(i).first) << std::endl;
            
            //std::cout << i << " " << ((int)alt.at(i) == (ve.at(i).first)) << std::endl;
            if ((int)alt.at(i) != ve.at(i).first) ++tot_mismatch;
        }
    }
    return tot_mismatch;
    
}


void print_alt(std::vector<size_t> vec) {
    for (auto i:vec) {
        std::cout << i << std::endl;
    }
}
    




int main(){
std::string astring {"ATGCGGGGGGGGGGTCAGACCCCCCTA"};


std::string st2 {"TATATATGCCGATCGGGATCCAT"};
std::string st3 {"TATATATACGATTCGATCGTAA"};
Timer ti;
ti.reset();
    
//Genome gen = Genome("/home/dan/pcr3/COH1.fasta");
ti.reset();
GenomeIndex gi2 {st2, 5};

//GenomeIndex gi2 {st2, 5};


//std::cout << ti.elapsed() << std::endl;
std::cout << "buckets" << gi2.get_index().bucket_count() << std::endl;
GenomeIndexSearch gs {st2, gi2};

std::string ss {"ATATGCCGA"};
gs.search_string(ss,0, 8);
gs.display_fwd_hits();

auto res = alt_find(st2, ss);
auto mism = compare_hits(res, gs);
std::cout << " Mismatches = "<<mism<<std::endl;

//std::cout << "Map size: " << estimate_memory_size(gi.get_index()) << std::endl;

return 0;
}
