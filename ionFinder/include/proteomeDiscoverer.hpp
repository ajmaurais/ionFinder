//
// proteomeDiscoverer.hpp
// ionFinder
// -----------------------------------------------------------------------------
// MIT License
// Copyright 2021 Aaron Maurais
// -----------------------------------------------------------------------------
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
// -----------------------------------------------------------------------------
//

#ifndef proteomeDiscoverer_hpp
#define proteomeDiscoverer_hpp

#include <fstream>
#include <vector>
#include <map>
#include <regex>
#include <algorithm>

#include <sqlite3.h>

#include <utils.hpp>

namespace pd{
    class Formula;
    class Peptide;
    class PSM;
    class AminoAcid;
    class Modification;

    class Formula{
    public:
        enum class ATOM {C, C13, H, D, Br, Cl, N, N15, O, O18, P, S, Se, Last, First = C};
        typedef std::map<ATOM, std::pair<double, double> > AtomMassType;

        static ATOM pd_str_to_atom(const std::string& s);
        static std::string atom_to_str(ATOM a);

    protected:
        std::map<ATOM, int> _formula;
        static AtomMassType _atomMasses;

    public:
        Formula() = default;

        void clear();
        void add(ATOM a, int count);
        void set(ATOM a, int count);
        void add(const Formula& rhs);
        void set(const Formula& rhs);

        double calcMass(bool mono = true) const;
        double calcAvgMass() const;
        double calcMonoMass() const;
        std::string calcFormula() const;
    };

    class AminoAcid: public Formula{
    private:
        std::string _name;
        std::string _oneLetterCode;
        std::string _threeLetterCode;
    private:
        double _monoMass;
        double _avgMass;
    public:
        AminoAcid() : Formula() {
            clear();
        }

        // modifers
        void clear();
        void setName(const std::string &name) {
            _name = name;
        }
        void setOneLetterCode(const std::string &oneLetterCode) {
            _oneLetterCode = oneLetterCode;
        }
        void setThreeLetterCode(const std::string &threeLetterCode) {
            _threeLetterCode = threeLetterCode;
        }
        void setAvgMass(double avgMass) {
            _avgMass = avgMass;
        }
        void setMonoMass(double monoMass) {
            _monoMass = monoMass;
        }

        // getters
        const std::string &getName() const {
            return _name;
        }
        const std::string &getOneLetterCode() const {
            return _oneLetterCode;
        }
        const std::string &getThreeLetterCode() const {
            return _threeLetterCode;
        }
        double getMonoMass() const {
            return _monoMass;
        }
        double getAvgMass() const {
            return _avgMass;
        }

    };

    inline Formula::ATOM operator++(Formula::ATOM& x ){
        return x = (Formula::ATOM)(((int)(x) + 1));
    }

    class Modification: public Formula{
    private:
        std::string _name;
        std::string _abbreviation;
    public:
        Modification(): Formula(){
            clear();
        }

        void clear();
        void setName(std::string s){
            _name = s;
        }
        void setAbbreviation(std::string s){
            _abbreviation = s;
        }

        // getters
        std::string getName() const{
            return _name;
        }
        std::string getAbbreviation() const{
            return _abbreviation;
        }
        std::string str() const;
    };

    class Residue{
    private:
        AminoAcid _aminoAcid;
        Modification _modification;
        bool _isModified;
    public:
        Residue(){
            _isModified = false;
        };
        explicit Residue(const AminoAcid& aa){
            _aminoAcid = aa;
            _isModified = false;
        }

        void setModification(const Modification& m){
            _isModified = true;
            _modification = m;
        }

        const AminoAcid& getAminoAcid() const{
            return _aminoAcid;
        }
        const Modification& getModification() const{
            return _modification;
        }
        AminoAcid& getAminoAcid(){
            return _aminoAcid;
        }
        Modification& getModification(){
            return _modification;
        }
        bool isModified() const{
            return _isModified;
        }

        std::string str() const;

    };

    class PeptideSequence: public Formula{
    private:
        std::vector<Residue> _residues;
        std::string sequence;

    public:
        PeptideSequence(): Formula(){}

        void clear();
        void setSequence(std::string, const std::map<std::string, AminoAcid>&);
        void addModification(const Modification&m, int index);

        std::string str() const;
        const std::string& getSequence() const{
            return sequence;
        }
//        const std::vector<Residue>& getResidues() const{
//            return _residues;
//        }
    };

    class Peptide{
    private:
        std::string _sequence;
        Formula _formula;
    public:

    };

    class PSM{
    private:
        std::string _formula;
        std::string _parentProtein;
        std::string _parentID;
        std::string _parentDescription;
        std::string _sampleName;

    };

    /**
    * Parse binary blob to a vector of numerical values.
    * PD also has an extra byte at the end of each value depicting whether the values is NULL.
    * \p isNA is populated with those flags.
    * @param arr vector of numerical values to populate.
    * @param isNA vector of NA flags.
    * @param data Data bits to parse.
    * @param dataLength Length of \p data in number of bytes.
    */
    template<class T>
    void getArray(std::vector<T>& arr, std::vector<bool>& isNA, unsigned char* data, size_t dataLength)
    {
        arr.clear();
        isNA.clear();

        T tmp;
        bool na;
        unsigned int valueSize = sizeof(T);
        unsigned int blockSize = valueSize + 1;
        size_t nValues = dataLength / blockSize;
        for(size_t i = 0; i < nValues; i++){
            memcpy(&tmp, data + i * blockSize, valueSize);
            memcpy(&na, data + i * blockSize + valueSize, 1);
            arr.push_back(tmp);
            isNA.push_back(!na);
        }
    }

    bool readPD(std::string fname, std::vector<Peptide>& peptides, std::vector<PSM>& PSMs);
    void parseAAFormula(const std::string& s, Formula& f);
    void parseModFormula(const std::string& s, Formula& f);
    void parseResidueModification(const std::string& s, Modification& m, const std::regex& re);
    std::string my_sqlite3_column_string(sqlite3_stmt* s, int i);
}

#endif
