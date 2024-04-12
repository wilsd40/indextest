#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>
#include "genome.h"
#include <climits>
#include <cctype>
#include <chrono>
#include <numeric>
#include <optional>
#include <utility>
#include <cassert>
#include <algorithm>
#include <map>
#include <stdexcept>




enum class Direction {FWD=0, COMP=1, REV_COMP=2};

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
    friend class GenomeIndexSearch;
    
    // note: ksize of 5 is a good compromise between speed of index creation (almost instantaneous) and search.
    // at least for a bacterial-sized search string.

using string_map = std::unordered_map<std::string, std::vector<int>>;
// this map provides constant access time, kmers are the string, and positions are in the vector.
using result_vector = std::vector<std::pair<int, int>>;
// results for query: [start position, length of identity].


private:
    string_map mp;
    string_map mp_comp;
    string_map mp_rev_comp;

    
public:
    int ksize;
    
GenomeIndex(const std::string& te1, int ksize1) {
    make_index(te1, ksize1, mp);
    make_index(te1, ksize1, mp_comp);
    make_index(te1, ksize1, mp_rev_comp);
    this->ksize = ksize;
    }
    
GenomeIndex(const Genome& ge, int ksize1){
    make_index(ge.get_genome(), ksize1, mp);
    make_index(ge.get_complement(), ksize1, mp_comp);
    make_index(ge.get_rev_comp(), ksize1, mp_rev_comp);
    this->ksize = ksize1;
}
    
void make_index(const std::string& te, int ksize, string_map& targ) {
    for (size_t i {0}; i < te.size() - ksize + 1; ++i) {
        std::string kmer = te.substr(i, ksize);
        targ[kmer].push_back(i); // no need the create an entry if one doesn't exist, it is already initialized
    }
    
}

// consider searching index here?


const string_map& get_index() const{
    return mp;
}

const string_map& get_comp_index() const {
    return mp_comp;
}
const string_map& get_rev_comp_index() const {
    return mp_rev_comp;
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
friend class Genome;
friend class GenomeIndex;
using string_map = std::unordered_map<std::string, std::vector<int>>;
// this map provides constant access time, kmers are the string, and positions are in the vector.
using result_vector = std::vector<std::pair<int, int>>;
// results for query: [start position, length of identity].    

private:
    const Genome& gen; // will be passed in via constructor
    const int ksize; // passed into constructor
    GenomeIndex gin; // will be created in constructor


    std::string quer;
    
    int min_match;
    int direction; // 0, fwd, 1, rev, 2, don't know
    
    result_vector fwd_result;
    result_vector rev_result;
    result_vector result_return;


    
public:
    static int validate_k_ge(const Genome& gea, int ks) {
        if ((int)gea.get_genome().size() < ks) {
            throw std::invalid_argument("ksize must less than genome size");
        }
        return ks;
    }
    static int validate_k_string(const std::string& gea, int ks) {
        if ((int)gea.size() < ks) {
            throw std::invalid_argument("ksize must less than genome size");
        }
        return ks;
    }
        
    std::map<char, char> comp = {{'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'}} ;
   

     
    GenomeIndexSearch(const Genome& gearg, int ks): gen {gearg}, ksize {validate_k_ge(gearg.get_genome(), ks)},
        gin {GenomeIndex(gearg.get_genome(), ksize)} {
        std::cout << "ksize = "<<ksize << std::endl;
        std::cout << "The genome size is "<<gen.get_genome().size() << std::endl;
    }
    
    //GenomeIndexSearch(std::string& gearg, const GenomeIndex& giarg, std::optional<int> arg3 = std::nullopt): ge {gearg}, gi {giarg.get_index()},
    //    ksize {giarg.get_ksize()}{
    //        std::cout << "ksize = "<<ksize << std::endl;
        // constructor for using string.
        // need to make complements.
     


    //}
// workhorse function for finding seed strings in index.  returns a vector of positions for which the string can be found.    
    std::vector<int> query_index(std::string query, Direction dir) {
        // make sure the check that query was length of index.
        std::string found;
        std::vector<int>rv;
        string_map& gioption {gin.mp}; // need to initialize it temporarily. not sure if this is the best
        // to do this but don't want to make a copy and need to declare it.
        
        
        switch (dir) {
            case Direction::FWD:
                // gioption already set to forward
                break;
            case Direction::COMP:
                gioption = gin.mp_comp;
                break;
            case Direction::REV_COMP:
                gioption = gin.mp_rev_comp;
                break;
        }
        
        auto it1 = gioption.find(query);
        if (it1 == gioption.end()) {
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
    
void sort_pairs(result_vector& vp) {
    std::sort(vp.begin(), vp.end(), [](auto val1, auto val2) {return ((val1.second != val2.second) ? (val1.second > val2.second) : (val1.first < val2.first));});    
}
void separate_pairs(const result_vector& vp, std::vector<int>& vec1, std::vector<int>& vec2){
    // this still need to be adapted (copy and pasted from other code; need to be adapted to class structure
    
    vec1.clear();
    vec2.clear();
    std::transform(vp.begin(), vp.end(), std::back_inserter(vec1), [] (const std::pair<int, int> p) {return p.first;});
    std::transform(vp.begin(), vp.end(), std::back_inserter(vec2), [] (const std::pair<int, int> p) {return p.second;});
}

// function to find primer in template.  calls query_index
bool search_string(std::string query, Direction dir, int min_match) {
    // this method will be called by another method to specify dir.  
    // search the subject with the query 
    // dir : 0 for forward primer (will start at 3' end and work backward on fwd strand
    //     : 1 for reverse primer (will make reverse complement, start at 5' end and work forward on fwd strand
    //int end_nuc {-1};
    
    //need to check min_match vs query.
    std::cout << (int)query.length() << ": ksize :" <<ksize << std::endl;
    assert ((int)query.length() >= ksize);
    assert (min_match >= ksize);
    if (dir == Direction::FWD) { // forward primer
        std::cout << "Searching in foward direction" << std::endl;
        std::string ge=gen.get_genome();
        //std::cout << "Genome is "  << ge << std::endl;
        fwd_result.clear();
        std::cout << ge.size() << " genome size" << std::endl;
        
        //std::cout << (int)query.length() << std::endl;
        int start_nuc {(int)query.length() - ksize}; // 
        
        
        std::cout << "Start nuc: "<<start_nuc<<std::endl;
 
        

        std::string qu_nucs = query.substr(start_nuc, ksize);
        std::cout << "Will start search with: "<< qu_nucs << std::endl;

        auto indexstart = query_index(qu_nucs, dir);
        if (indexstart.size() == 1 && indexstart.at(0) == -1) {
            std::cout << "search string not found" << std::endl;
            fwd_result.push_back(std::pair<int, int>{-1, -1});
            return 1;
            
        }
        std::cout << "Index size "<<indexstart.size() << std::endl;
        for (auto fi:indexstart) {
        //auto fi = starttry.at(0);
        
        //for (auto fi:query_res){
            
            
            std::cout << "Starting at position "<<fi << std::endl;
           
        
   
            std::string result_string;
            int xtend {0};
            std::cout << "conditions" << (start_nuc - xtend - 1) << " " << (fi - xtend - 1) << std::endl;
            // shouldn't have to check 3' bounds
            //bool found = false;
            while (((start_nuc - xtend - 1) >= 0 ) && (fi - xtend - 1) >=0){
                std::cout << "------" << std::endl;
                //std::cout << "before: bounds checks : start_nuc-xtend-1: " << start_nuc-xtend-1 << " fi - xtend - 1:  " << (fi - xtend - 1)<<std::endl;
                //std::cout << xtend << " comparing " << query.at(start_nuc - xtend - 1) << " to " << ge.at(fi - xtend - 1) << std::endl;
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
    sort_pairs(fwd_result);
    
    
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
    const std::string& ge {gen.genome};
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

void ToUpper(std::string& st) {
    std::transform(st.begin(), st.end(), st.begin(), ::toupper);
}


    




int main(){
    
// AAGCTACCTGCTAGGGCGAC
    
Genome g2 = Genome("/home/dan/pcr3/test.fasta", false);
std::cout << g2.get_size() << std::endl;
std::cout << g2.get_genome() << std::endl;
GenomeIndexSearch gs {g2, 7};
auto res = gs.query_index("AGGGCGA", Direction::FWD);
for (auto i:res) {
    std::cout << i << " ";
}

gs.search_string("TAGGGCGA", Direction::FWD, 7);
gs.display_fwd_hits();



//std::string astring {"ATGCGGGGGGGGGGTCAGACCCCCCTA"};

//GenomeIndexSearch gs {g2, 7};

/*
std::string st2 {"TATATATGCCGATCGGGATCCAT"};
std::string st3 {"TATATATACGATTCGATCGTAA"};
Timer ti;
ti.reset();
    
Genome gen = Genome("/home/dan/pcr3/COH1.fasta");
Genome g2 = Genome("/home/dan/pcr3/test.fasta");
ti.reset();
std::cout << g2.get_genome() << std::endl;
std::cout << g2.get_complement() << std::endl;
std::cout << g2.get_rev_comp() << std::endl;*/
//GenomeIndex gi2 {gen, 5};

//GenomeIndex gi2 {st2, 5};


//std::cout << ti.elapsed() << std::endl;
//std::cout << "buckets" << gi2.get_index().bucket_count() << std::endl;
//GenomeIndexSearch gs {gen, gi2};

//std::string ss {"ATATGCCGA"};
//gs.search_string(ss,0, 7);
//gs.display_fwd_hits();

//auto res = alt_find(st2, ss);
//auto mism = compare_hits(res, gs);
//std::cout << " Mismatches = "<<mism<<std::endl;

//std::cout << "Map size: " << estimate_memory_size(gi.get_index()) << std::endl;

return 0;
}
