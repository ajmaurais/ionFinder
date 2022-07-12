//
// params.cpp
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

#include <ionFinder/params.hpp>

/**
 Calculates number of supported threads using std::thread::hardware_concurrency(). <br><br>
 By default returns result of thread::hardware_concurrency() / 2.<br>
 If thread::hardware_concurrency returns 0, indicating an error, 1 thread is used.
 \return number of supported threads
 */
unsigned int IonFinder::Params::computeThreads()
{
    unsigned int ret = std::thread::hardware_concurrency();
    if(ret == 0){
        std::cerr << "\nError detecting hardware_concurrency. Only 1 thread being used." << NEW_LINE;
        return 1;
    }
    return ret;
}

IonFinder::Params::InputFileType IonFinder::Params::strToInputFileType(std::string s) {
    if(s == IonFinder::DTAFILTER_INPUT_STR) return InputFileType::DTAFILTER;
    if(s == IonFinder::TSV_INPUT_STR) return InputFileType::TSV;
    if(s == IonFinder::MZ_IDENT_ML_STR) return InputFileType::MZ_IDENT_ML;
    throw std::invalid_argument("'" + s + "' is not a valid InputFileType!");
}

std::string IonFinder::Params::inputFileTypeToStr(IonFinder::Params::InputFileType it) {
    if(it == InputFileType::DTAFILTER) return IonFinder::DTAFILTER_INPUT_STR;
    if(it == InputFileType::TSV) return IonFinder::TSV_INPUT_STR;
    if(it == InputFileType::MZ_IDENT_ML) return IonFinder::MZ_IDENT_ML_STR;
    throw std::runtime_error("No string conversion for InputFileType!");
}

/**
 Parses command line arguments and stores in Params object
 \pre current working directory exists
 \param argc argc from main
 \param argv argv from main
 \return true if args were successfully read.
 */
bool IonFinder::Params::getArgs(int argc, const char* const argv[])
{
    //Should program continue if no filter files are found in a dir?
    bool force = false;
    _wd = utils::pwd();
    assert(utils::dirExists(_wd));

    for(int i = 1; i < argc; i++)
    {
        if(!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help"))
        {
            displayHelp();
            return false;
        }
        if(!strcmp(argv[i], "-d") || !strcmp(argv[i], "--dir"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            _wd = utils::absPath(argv[i]);
            wdSpecified = true;
            if(!utils::dirExists(_wd))
            {
                std::cerr << "Specified directory does not exist." << NEW_LINE;
                return false;
            }
            continue;
        }
        if(!strcmp(argv[i], "-i") || !strcmp(argv[i], "--inputMode"))
        {
            if(!utils::isArg(argv[++i])) {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            try {
                _inputMode = strToInputFileType(argv[i]);
            }
            catch (const std::invalid_argument& e) {
                std::cerr << e.what() << NEW_LINE;
                return false;
            }
            continue;
        }
        if(!strcmp(argv[i], "-o") || !strcmp(argv[i], "--ofname"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            ofname = argv[i];
            continue;
        }
        if(!strcmp(argv[i], "-dta"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            _dtaFilterBase = argv[i];
            continue;
        }
        if(!strcmp(argv[i], "-rev"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            if(!(!strcmp(argv[i], "0") || !strcmp(argv[i], "1")))
            {
                std::cerr << argv[i] << base::PARAM_ERROR_MESSAGE << argv[i-1] << std::endl;
                return false;
            }
            _includeReverse = std::stoi(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "-m") || !strcmp(argv[i], "--modFilter"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            if(!(!strcmp(argv[i], "0") || !strcmp(argv[i], "1") || !strcmp(argv[i], "2")))
            {
                std::cerr << argv[i] << base::PARAM_ERROR_MESSAGE << argv[i-1] << std::endl;
                return false;
            }
            _modFilter = inputFiles::intToModFilter(std::stoi(argv[i]));
            continue;
        }
        if(!strcmp(argv[i], "--fastaFile"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            _fastaFile = utils::absPath(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "-u") || !strcmp(argv[i], "--peptideUID"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            if(!(!strcmp(argv[i], "0") || !strcmp(argv[i], "1")))
            {
                std::cerr << argv[i] << base::PARAM_ERROR_MESSAGE << argv[i-1] << std::endl;
                return false;
            }
            _printPeptideUID = std::stoi(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "-I") || !strcmp(argv[i], "--printInt"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            if(!(!strcmp(argv[i], "0") || !strcmp(argv[i], "1")))
            {
                std::cerr << argv[i] << base::PARAM_ERROR_MESSAGE << argv[i-1] << std::endl;
                return false;
            }
            _printIonIntensity = std::stoi(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "-p") || !strcmp(argv[i], "--printSpectra"))
        {
            _printSpectraFiles = true;
            continue;
        }
        if(!strcmp(argv[i], "--calcNL"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            if(!(!strcmp(argv[i], "0") || !strcmp(argv[i], "1")))
            {
                std::cerr << argv[i] << base::PARAM_ERROR_MESSAGE << argv[i-1] << NEW_LINE;
                return false;
            }
            _calcNL = std::stoi(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "-l") || !strcmp(argv[i], "--lossMass"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            //TODO: Maybe add option to search for multiple losses
            _neutralLossMass = std::stod(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "--modMass"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            _modMass = std::stod(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "-g") || !strcmp(argv[i], "--groupMod"))
        {
            if(!utils::isArg(argv[++i])) {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            if(!(!strcmp(argv[i], "0") || !strcmp(argv[i], "1")))
            {
                std::cerr << argv[i] << base::PARAM_ERROR_MESSAGE << argv[i-1] << NEW_LINE;
                return false;
            }
            _groupMod = std::stoi(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "-mt") || !strcmp(argv[i], "--matchTolerance"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            matchTolerance = std::stod(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "--matchType"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            _matchType = strToMatchType(std::string(argv[i]));
            if(_matchType == MatchType::UNKNOWN){
                std::cerr << argv[i] << " is an unknown MatchType!" << NEW_LINE;
                return false;
            }
            continue;
        }
        if(!strcmp(argv[i], "--citStats"))
        {
            _modMass = ionFinder::CIT_MOD_MASS;
            _neutralLossMass = ionFinder::CIT_NL_MASS;
            _ambigiousResidues = ionFinder::CIT_AMB_RESIDUES;
            _includeCTermMod = ionFinder::CIT_INCLUDE_C_TERM_MOD;
            _calcNL = ionFinder::CIT_CALC_NL;
            ofname = ionFinder::PEPTIDE_CIT_STATS_OFNAME;
            continue;
        }
        if(!strcmp(argv[i], "--cTermMod"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            _includeCTermMod = std::stoi(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "--isoAA"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            _ambigiousResidues = std::string(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "-minC"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            minFragCharge = std::stoi(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "-maxC"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            maxFragCharge = std::stoi(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "-minMZ"))
        {
            if(!utils::isArg(argv[++i])){
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            minMZ = std::stod(argv[i]);
            minMzSpecified = true;
            continue;
        }
        if(!strcmp(argv[i], "-maxMZ"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            maxMZ = std::stod(argv[i]);
            maxMzSpecified = true;
            continue;
        }
        if(!strcmp(argv[i], "-minLabInt"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            minLabelIntensity = std::stoi(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "-minNlLabInt"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            _minNlLabelIntensity = std::stod(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "-n") || !strcmp(argv[i], "--artifactNLIntPerc"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            _artifactNLIntFrac = std::stod(argv[i]) / 100; // Convert from percentage to fraction here.
            continue;
        }
        if(!strcmp(argv[i], "--labelArtifactNL"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            if(!(!strcmp(argv[i], "0") || !strcmp(argv[i], "1")))
            {
                std::cerr << argv[i] << base::PARAM_ERROR_MESSAGE << argv[i-1] << NEW_LINE;
                return false;
            }
            _labelArtifactNL = std::stoi(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "-minInt"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            minIntensity = std::stod(argv[i]);
            minIntensitySpecified = true;
            continue;
        }
        if(!strcmp(argv[i], "-minSNR"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            _minSNR = std::stod(argv[i]);
            minSNRSpecified = true;
            continue;
        }
        if(!strcmp(argv[i], "--snrConf"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            _snrConf = std::stod(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "-y") || !strcmp(argv[i], "--height"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            plotHeight = std::stod(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "-w") || !strcmp(argv[i], "--width"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            plotWidth = std::stod(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "-mmComp"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            if(!(!strcmp(argv[i], "intensity") || !strcmp(argv[i], "mz")))
            {
                std::cerr << argv[i] << base::PARAM_ERROR_MESSAGE << argv[i-1] << NEW_LINE;
                return false;
            }
            multipleMatchCompare = std::string(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "--incAllIons"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            if(!(!strcmp(argv[i], "0") || !strcmp(argv[i], "1")))
            {
                std::cerr << argv[i] << base::PARAM_ERROR_MESSAGE << argv[i-1] << NEW_LINE;
                return false;
            }
            includeAllIons = std::stoi(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "--printSmod"))
        {
            if(!writeSmod(_wd))
                std::cerr << "Could not write new smod file!" << NEW_LINE;
            return false;
        }
        if(!strcmp(argv[i], "--smod"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage();
                return false;
            }
            _smodFile = utils::absPath(argv[i]);
            smodSpecified = true;
            continue;
        }
        if(!strcmp(argv[i], "--nThread"))
        {
            if(!utils::isArg(argv[++i]))
            {
                usage(IonFinder::ARG_REQUIRED_STR + argv[i-1]);
                return false;
            }
            _numThread = std::stoi(argv[i]);
            continue;
        }
        if(!strcmp(argv[i], "--parallel"))
        {
            _numThread = computeThreads();
            continue;
        }
        if(!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbose"))
        {
            verbose = true;
            continue;
        }
        if(!strcmp(argv[i], "--version"))
        {
            printVersion(std::cout);
            printGitVersion(std::cout);
            return false;
        }
        if(!strcmp(argv[i], "-f")){
            force = true;
        }
        else if(utils::isFlag(argv[i])){
            std::cerr << argv[i] << " is an invalid argument." << NEW_LINE;
            usage();
            return false;
        }
        else{ //we are in input dirs
            while(i < argc){
                if(utils::isFlag(argv[i])){
                    usage();
                    return false;
                }
                _inDirs.emplace_back(argv[i++]);
                _inDirSpecified = true;
            }
        }
        //end of else
    }//end of for i

    //fix options
    if(_wd[_wd.length() - 1] != '/')
        _wd += "/";
    if(_inputMode == InputFileType::DTAFILTER){
        if(!getFlist(force)){
            std::cerr << "Could not find DTAFilter-files!" << NEW_LINE;
            return false;
        }
    }
    else if(_inputMode == InputFileType::TSV || _inputMode == InputFileType::MZ_IDENT_ML) {
        if(_inDirs.empty()) {
            std::cerr << "ERROR: Input file name is required when using " + inputFileTypeToStr(_inputMode) + " input mode!\n";
            usage();
            return false;
        }
    }

    return true;
} // end of getArgs

/**
 Searches all directories in _inDirs for DTAFilter files.
 If _inDirs is empty, current working directory is used.
 \return true if > 1 filter file was found, else false
 */
bool IonFinder::Params::getFlist(bool force)
{
    if(_inDirs.empty()){
        _inDirs.push_back(_wd);
        _wd = utils::parentDir(_wd);
    }
    for(auto& _inDir : _inDirs)
    {
        std::string fname = (_inDirSpecified ? (_wd + _inDir) : _inDir) + ("/" + _dtaFilterBase);
        if(utils::fileExists(fname)){
            _filterFiles[utils::baseName(_inDir)] = fname;
        }
        else{
            if(force) std::cerr << "WARN: ";
            else std::cerr << "ERROR: ";
            std::cerr << "No filter file found in: " << _inDir << NEW_LINE;
            if(force) return false;
        }
    }
    return !_filterFiles.empty();
}

/**
 * Print program version and git repo info to stream.
 * \param out Stream to print to.
 */
void IonFinder::Params::printVersion(std::ostream& out) {
    out << "ionFinder v" << PROG_VERSION_MAJOR << "."
        << PROG_VERSION_MINOR << "."
        << PROG_VERSION_PATCH << NEW_LINE;
}
