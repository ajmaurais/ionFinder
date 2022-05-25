//
// proteomeDiscoverer.cpp
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

#include <proteomeDiscoverer.hpp>

namespace pd {

    Formula::AtomMassType Formula::_atomMasses{
        {Formula::ATOM::C, {12.011, 12} },
        {Formula::ATOM::H, {1.008, 1.00783}},
        {Formula::ATOM::O, {15.999, 15.99491}},
        {Formula::ATOM::O18, {17.99916, 17.99916}},
        {Formula::ATOM::N, {14.007, 14.00307}},
        {Formula::ATOM::S, {32.06, 31.97207}},
        {Formula::ATOM::P, {30.97376, 30.97376}},
        {Formula::ATOM::N15, {15.00011, 15.00011}},
        {Formula::ATOM::D, {2.0141, 2.0141}},
        {Formula::ATOM::C13, {13.00335, 13.00335}},
        {Formula::ATOM::Se, {78.96, 79.91652}},
        {Formula::ATOM::Cl, {35.45, 34.96885}},
        {Formula::ATOM::Br, {79.904, 78.91834}}
    };

    //! Convert PD atom string to ATOM
    Formula::ATOM Formula::pd_str_to_atom(const std::string &s) {
        if(s == "C") return Formula::ATOM::C;
        else if(s == "H") return Formula::ATOM::H;
        else if(s == "O") return Formula::ATOM::O;
        else if(s == "18O") return Formula::ATOM::O18;
        else if(s == "N") return Formula::ATOM::N;
        else if(s == "S") return Formula::ATOM::S;
        else if(s == "P") return Formula::ATOM::P;
        else if(s == "15N") return Formula::ATOM::N15;
        else if(s == "2H") return Formula::ATOM::D;
        else if(s == "13C") return Formula::ATOM::C13;
        else if(s == "Se") return Formula::ATOM::Se;
        else if(s == "Cl") return Formula::ATOM::Cl;
        else if(s == "Br") return Formula::ATOM::Br;
        else throw std::runtime_error(s + " is an unknown atom!");
    }

    //! Convert ATOM to string in ionFinder formula format
    std::string Formula::atom_to_str(Formula::ATOM a) {
        switch(a) {
            case Formula::ATOM::C: return "C"; break;
            case Formula::ATOM::C13: return "(13)C"; break;
            case Formula::ATOM::H: return "H"; break;
            case Formula::ATOM::D: return "D"; break;
            case Formula::ATOM::Br: return "Br"; break;
            case Formula::ATOM::Cl: return "Cl"; break;
            case Formula::ATOM::N: return "N"; break;
            case Formula::ATOM::N15: return "(15)N"; break;
            case Formula::ATOM::O: return "O"; break;
            case Formula::ATOM::O18: return "(18)O"; break;
            case Formula::ATOM::P: return "P"; break;
            case Formula::ATOM::S: return "S"; break;
            case Formula::ATOM::Se: return "Se"; break;
            default: throw std::runtime_error("Unknown atom!");
        }
    }

    void Formula::clear() {
        _formula.clear();
    }

    void Modification::clear(){
        _name = "";
        _abbreviation = "";
    }

    //

    void Formula::add(Formula::ATOM a, int count) {
        if(_formula.find(a) == _formula.end())
            set(a, count);
        else _formula[a] += count;
    }

    void Formula::set(Formula::ATOM a, int count) {
        _formula[a] = count;
    }

    void Formula::set(const Formula& rhs){
        clear();
        for(const auto& it: rhs._formula)
            set(it.first, it.second);
    }

    void Formula::add(const Formula& rhs){
        for(const auto& it: rhs._formula)
            add(it.first, it.second);
    }

    /**
     * Use Formula to calculate mass.
     * @param mono Should the mono mass be calculated? If false, the avg mass is calculated instead.
     * @return mass
     */
    double Formula::calcMass(bool mono) const{
        double ret = 0;
        for(const auto& it: _formula)
            ret += (mono ? _atomMasses[it.first].second : _atomMasses[it.first].first) * it.second;
        return ret;
    }

    //! Use Formula to calculate avg mass
    double Formula::calcAvgMass() const {
        return calcMass(false);
    }

    //! Use Formula to calculate mono mass
    double Formula::calcMonoMass() const {
        return calcMass(true);
    }

    //! Get string formatted formula
    std::string Formula::calcFormula() const {
        std::string ret = "";
        for(auto a = ATOM::First; a != ATOM::Last; ++a){
            if(_formula.find(a) != _formula.end()) {
                if(_formula.at(a) > 0) {
                    ret += atom_to_str(a);
                    if (_formula.at(a) > 1)
                        ret += std::to_string(_formula.at(a));
                }
            }
        }
        return ret;
    }

    void AminoAcid::clear() {
        Formula::clear();
        _name = "";
        _oneLetterCode = "";
        _threeLetterCode = "";
        _monoMass = 0;
        _avgMass = 0;
    }

    /**
     * Parse formula given in FoundModificationsTargetPsms table.
     * @param s Raw formula string.
     * @param f Empty Formula to populate.
     */
    void parseModFormula(const std::string& s, Formula& f)
    {
        std::vector<std::string> elems;
        utils::split(s, ' ', elems);
        f.clear();

        std::regex re(R"(([A-Z0-9]+)(?:\(([0-9]+)\))?)");
        std::smatch matches;
        for(const auto& e: elems){
            if(!std::regex_search(e, matches, re))
                throw std::runtime_error("Could not parse atom element: " + e);
            f.add(Formula::pd_str_to_atom(matches[1]), matches[2] == "" ? 1 : std::stoi(matches[2]));
        }
    }

    /**
     * Parse formula given in AminoAcids table.
     * @param s Raw formula string.
     * @param f Empty Formula to populate.
     */
    void parseAAFormula(const std::string& s, Formula& f)
    {
        std::vector<bool> indices(s.length(), false); // keep track of indices which have been covered by regex

        std::regex re(R"(([A-Z][a-z]?)([\d]+)?)");
        for(auto it = std::sregex_iterator(s.begin(), s.end(), re); it != std::sregex_iterator(); ++it){
            f.add(Formula::pd_str_to_atom((*it)[1]), (*it)[2] == "" ? 1 : std::stoi((*it)[2]));

            // Record regex traversed indices
            unsigned int mBegin = it->position();
            unsigned int mEnd = mBegin + it->length();
            for(unsigned int i = mBegin; i < mEnd; i++)
                indices[i] = true;
        }

        // check that all string indices were traversed
        const bool all_traversed = std::all_of(indices.begin(), indices.end(), [](const bool x){ return x; });
        if(!all_traversed) {
            std::cerr << NEW_LINE << "Regex error parsing formula:" << NEW_LINE << "\t";
            for(const auto& c: s){
                std::cerr << c << " ";
            }
            std::cerr << NEW_LINE << "\t";
            for(const auto& b: indices){
                std::cerr << int(b) << " ";
            }
            std::cerr << NEW_LINE;
            throw std::runtime_error("Unable to parse formula: " + s);
        }
    }

    std::string my_sqlite3_column_string(sqlite3_stmt *s, int i) {
        return std::string(reinterpret_cast<const char*>(sqlite3_column_text(s, i)));
    }

    std::string Modification::str() const{
        std::string ret = "(";

        return ret;
    }

    std::string Residue::str() const {
        return _aminoAcid.getOneLetterCode() + _modification.str();
    }

    std::string PeptideSequence::str() const {
        std::string ret = "";

        return ret;
    }

    void PeptideSequence::clear() {
        Formula::clear();
        _residues.clear();
        sequence.clear();
    }

    /**
     * Populate sequence with Residue objects.
     * @param s string representation of peptide.
     * @param aa A map of amino acids objects used to build sequence.
     */
    void PeptideSequence::setSequence(std::string s, const std::map<std::string, AminoAcid>& aa) {
        clear();
        sequence = s;

        // update residues
        _residues.emplace_back(aa.at("Nte")); // Add n term
        for(auto c: s) _residues.emplace_back(aa.at(std::string(1, c))); // add residues
        _residues.emplace_back(aa.at("Cte")); // Add n term

        // Update formula
        // Skip modifications because that should be empty at this point
        for(auto r: _residues) {
            add(r.getAminoAcid());
        }
    }

    /**
     * Add Modification \p m at index \p index.
     * @param m Modification to add.
     * @param index Index of modification. The first amino acid is index 1. The n terminus is index 0.
     * The C terminus is index sequence.length() + 1.
     */
    void PeptideSequence::addModification(const Modification &m, int index) {
        _residues.at(index).setModification(m);
        add(m);
    }

    bool readPD(std::string fname, std::vector<Peptide> &peptides, std::vector<PSM> &PSMs)
    {
        //make sure PD file exists
        if (!utils::fileExists(fname)) {
            std::cerr << fname << " does not exist!" << NEW_LINE;
            return false;
        }

        // open database connection
        sqlite3* connection;
        sqlite3_open(fname.c_str(), &connection);
        sqlite3_stmt* stmt;
        int rc;

        // get amino acids
        std::map<std::string, AminoAcid> aminoAcids;
        std::string aa_query = R"(
            SELECT
                Name,               -- 0
                ThreeLetterCode,    -- 1
                OneLetterCode,      -- 2
                MonoisotopicMass,   -- 3
                AverageMass,        -- 4
                SumFormula          -- 5
            FROM AminoAcids
            WHERE SumFormula IS NOT "";
        )";
        sqlite3_prepare_v2(connection, aa_query.c_str(), -1, &stmt, nullptr);
        rc = sqlite3_step(stmt);
        if(rc != SQLITE_ROW){
            std::cerr << NEW_LINE << "ERROR: Could not retrieve list of amino acids!" << NEW_LINE;
            return false;
        }
        for(; rc == SQLITE_ROW; rc = sqlite3_step(stmt)) {
            AminoAcid temp;
            std::string oneLetter = my_sqlite3_column_string(stmt, 2);
            std::string threeLetter = my_sqlite3_column_string(stmt, 1);
            std::string key;
            if(threeLetter == "Nte" || threeLetter == "Cte") key = threeLetter; // use three letter for termini
            else key = oneLetter; // use one letter for residues
            temp.setThreeLetterCode(threeLetter);
            temp.setName(my_sqlite3_column_string(stmt, 0));
            temp.setOneLetterCode(my_sqlite3_column_string(stmt, 2));
            temp.setMonoMass(sqlite3_column_double(stmt, 3));
            temp.setAvgMass(sqlite3_column_double(stmt, 4));
            parseAAFormula(my_sqlite3_column_string(stmt, 5), temp);

            aminoAcids[key] = temp;
        }

        // get modifications
        std::map<std::string, Modification> modifications;
        std::string mod_query = R"(
            SELECT
                Abbreviation,   -- 0
                Name,           -- 1
                Substitution    -- 2
            FROM FoundModifications
            WHERE ModificationID IN
                (SELECT DISTINCT
                    FoundModificationsModificationID
                FROM FoundModificationsTargetPsms);
        )";
        sqlite3_prepare_v2(connection, mod_query.c_str(), -1, &stmt, nullptr);
        rc = sqlite3_step(stmt);
        if(rc != SQLITE_ROW){
            std::cerr << NEW_LINE << "ERROR: Could not retrieve modifications!" << NEW_LINE;
            return false;
        }
        for(; rc == SQLITE_ROW; rc = sqlite3_step(stmt)) {
            Modification temp;
            std::string s = my_sqlite3_column_string(stmt, 0);
            temp.setAbbreviation(s);
            std::string f = my_sqlite3_column_string(stmt, 2);
            parseModFormula(f, temp);
            temp.setName(my_sqlite3_column_string(stmt, 1));
            modifications[s] = temp;
        }

        // std::vector<double> diffs;
        // for(auto aa: aminoAcids){
        //     diffs.push_back(abs(aa.second.calcMonoMass() - aa.second.getMonoMass()));
        //     std::cout << "Mono calcd: " << aa.second.calcMonoMass() << ", Mono obsd:" << aa.second.getMonoMass() << NEW_LINE;
        //     //std::cout << "Mono calcd: " << aa.second.calcMonoMass() << ", Mono obsd:" << aa.second.getMonoMass() << NEW_LINE;
        // }
        // std::cout << "\nMax difference is: " << *std::max_element(diffs.begin(), diffs.end()) << NEW_LINE;

        // get PSMs
        std::string psm_query = R"(
            SELECT
                ids.PSMPeptideID,                   -- 0
                ids.PeptidePeptideID,               -- 1
                FirstScan AS ScanNum,               -- 2
                Sequence,                           -- 3
                Modifications,                      -- 4
                SpectrumFileName,                   -- 5
                StudyFileId,                        -- 6
                protIDs.ProteinID,                  -- 7
                protIDs.ProteinDescription,         -- 8
                protIDs.FullSequence,               -- 9
                protIDs.PositionsinMasterProteins,  -- 10
                tp.Charge,                          -- 11
                tp.XCorr,                           -- 12
                MassOverCharge,                     -- 13
                MasterScanNumbers AS PrecursorScan  -- 14
            FROM TargetPsms tp
            LEFT JOIN
                (SELECT
                 TargetPsmsPeptideID AS PSMPeptideID,
                 TargetPeptideGroupsPeptideGroupID AS PeptidePeptideID
            FROM TargetPsmsTargetPeptideGroups) AS ids
            ON ids.PSMPeptideID == tp.PeptideID
            LEFT JOIN
                (SELECT
                    MasterProteinAccessions AS ProteinID,
                    MasterProteinDescriptions AS ProteinDescription,
                    AnnotatedSequence AS FullSequence,
                    PositionsinMasterProteins,
                    PeptideGroupID AS PeptidePeptideID
                FROM TargetPeptideGroups) AS protIds
            ON ids.PeptidePeptideID == protIDs.PeptidePeptideID
            WHERE protIDs.ProteinID IS NOT NULL;
        )";
        sqlite3_prepare_v2(connection, psm_query.c_str(), -1, &stmt, nullptr);
        rc = sqlite3_step(stmt);
        if(rc != SQLITE_ROW){
            std::cerr << NEW_LINE << "ERROR: Could not retrieve modifications!" << NEW_LINE;
            return false;
        }
        std::regex modification_regex = std::regex(R"((?:([A-Z])([0-9]+)|([NC]-Term)(?:\([A-Za-z0-9]+\))?)\(([A-Za-z0-9+\-_ ]+)\))");
        PeptideSequence pepSeq;
        for(; rc == SQLITE_ROW; rc = sqlite3_step(stmt)) {
            pepSeq.clear();

            // parse sequence
            std::string seq = my_sqlite3_column_string(stmt, 3);
            pepSeq.setSequence(seq, aminoAcids);

            //parse modifications
            std::string mods = my_sqlite3_column_string(stmt, 4);
            if(!mods.empty()) {
                std::vector<std::string> elems;
                utils::split(mods, ';', elems);
                utils::trimAll(elems);
                for (const auto& m: elems) {
                    Modification mod;
                    std::smatch matches;
                    if (!std::regex_search(m, matches, modification_regex))
                        throw std::runtime_error("Could not parse modification: \"" + m + "\"\n");

                    // get modification object corresponding to name
                    std::string key = matches[4];
                    const auto mod_temp = modifications.find(key);
                    if(mod_temp == modifications.end())
                        throw std::runtime_error("Unknown modification: \"" + m + "\"\n");
                    mod = mod_temp->second;

                    // process differently depending on whether mod is on peptide terminus or residue
                    if(matches[3] == "") {  // mod is on residue
                        int number = std::stoi(matches[2]);
                        char res = std::string(matches[1])[0]; // This is ok because the regex only allows 1 char
                        if(pepSeq.getSequence()[number - 1] != res)
                            throw std::runtime_error("Modifications do not match sequence: " + seq + ", " + m);
                        pepSeq.addModification(mod, number);
                    }
                    else{ // mod is on terminus
                        std::string terminus = matches[3];
                        int number;
                        if(terminus == "N-term") number = 0;
                        else if(terminus == "C-term") number = int(pepSeq.getSequence().length()) + 1;
                        else throw std::runtime_error("Unknown peptide terminus: \"" + terminus + "\"\n");
                        pepSeq.addModification(mod, number);
                    }

                } // end for each modification
            } // end if !mods.empty()
        } //end for each SQLite row

        // get peptides

        // close DB connection
        sqlite3_finalize(stmt);
        sqlite3_close(connection);
        return true;
    }

//    void parseResidueModification(const std::string &s, Modification &m, const std::regex& re) {
//        std::smatch matches;
//        for(const auto& e: elems){
//            if(!std::regex_search(e, matches, re))
//                throw std::runtime_error("Could not parse atom element: " + e);
//    }


} // end of namespace
