//
// inputFiles.hpp
// ionFinder
// -----------------------------------------------------------------------------
// MIT License
// Copyright 2020 Aaron Maurais
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

#ifndef inputFiles_hpp
#define inputFiles_hpp

#include <vector>

#include <scanData.hpp>
#include <input_file_constants.hpp>

namespace inputFiles {

    class Scan;
    class InputFile;

    class Scan : public scanData::Scan{
    public:
        enum class MatchDirection{FORWARD, REVERSE};

    private:
        static std::string const REVERSE_MATCH;
        std::string _formula;
        std::string _parentProtein;
        std::string _parentID;
        std::string _parentDescription;
        MatchDirection _matchDirection;
        std::string _sampleName;
        bool _unique;

    public:
        Scan() : scanData::Scan(){
            _formula = "";
            _parentProtein = "";
            _parentID = "";
            _parentDescription = "";
            _unique = false;
            _matchDirection = MatchDirection::REVERSE;
        }

        //modifiers
        Scan& operator = (const Scan&);
        void setFormula(std::string s){
            _formula = s;
        }
        void setParentProtein(std::string s){
            _parentProtein = s;
        }
        void setParentID(std::string s){
            _parentID = s;
        }
        void setMatchDirection(MatchDirection m){
            _matchDirection = m;
        }
        void setSampleName(std::string s){
            _sampleName = s;
        }
        void setParentDescription(std::string s){
            _parentDescription = s;
        }
        void setUnique(bool boo){
            _unique = boo;
        }

        //properties
        std::string getFormula() const{
            return _formula;
        }
        std::string getParentProtein() const{
            return _parentProtein;
        }
        std::string getParentID() const{
            return _parentID;
        }
        MatchDirection getMatchDirection() const{
            return _matchDirection;
        }
        std::string getSampleName() const{
            return _sampleName;
        }
        std::string getParentDescription() const{
            return _parentDescription;
        }
        bool getUnique() const{
            return _unique;
        }
        static Scan::MatchDirection strToMatchDirection(std::string str);
    };

    class InputFile {
    public:
        enum class InputFileType {DTA_FILTER, TSV, MZ_IDENT_ML};
        //! Represents which peptides to read from input files
        enum class ModFilter { ONLY_MODIFIED, ALL, EXCLUDE_MODIFIED };
        static ModFilter intToModFilter(int);

        static InputFileType strToInputFileType(std::string);
        static std::string inputFileTypeToStr(InputFileType);

    protected:
        //!how will peptides to be searched for be supplied?
        InputFileType _fileType;
        //!Contains all scan input files to be read
        std::vector<std::pair<std::string, std::string> > _inputFiles;
        //! should reverse peptide matches be considered
        bool _includeReverse;
        //! Which modification statuses should be included in output?
        ModFilter _modFilter;
        //! The input file extension
        std::string _fileExtension;
        //! The input file basename
        std::string _fileBasename;

        virtual bool read(const std::string&, std::vector<Scan>&, const std::string&) const = 0;

    public:
        InputFile() {
            _fileType = InputFileType::DTA_FILTER;
            _includeReverse = false;
            _modFilter = ModFilter::ALL;
            _fileExtension = "";
            _fileBasename = "";
        }

        void setInputFileType(InputFileType ift){
            _fileType = ift;
        }
        void setIncludeReverse(bool includeReverse) {
            _includeReverse = includeReverse;
        }
        void setModFilter(ModFilter modFilter) {
            _modFilter = modFilter;
        }
        void setFileExtension(std::string e) {
            _fileExtension = e;
        }
        void setFileBasename(std::string n) {
            _fileBasename = n;
        }

        InputFileType getInputFileType() const{
            return _fileType;
        }
        bool getIncludeReverse() const {
            return _includeReverse;
        }
        ModFilter getModFilter() const {
            return _modFilter;
        }
        std::string getFileExtension() const {
            return _fileExtension;
        }
        std::string getFileBasename() const {
            return _fileBasename;
        }

        virtual bool findInputFiles(const std::vector<std::string>& inputArgs, std::string& wd) = 0;
        bool read(std::vector<Scan>&, bool verbose = false) const;
    };


}

#endif /* inputFiles_hpp */