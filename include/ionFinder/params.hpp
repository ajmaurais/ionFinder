//
//  params.hpp
//  citFinder
//
//  Created by Aaron Maurais on 12/9/18.
//  Copyright © 2018 Aaron Maurais. All rights reserved.
//

#ifndef params_hpp
#define params_hpp

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <thread>

#include <ionFinder/ionFinder.hpp>
#include <paramsBase.hpp>
#include <utils.hpp>

namespace IonFinder{
	
	//program file locations
	std::string const PROG_USAGE_FNAME = base::PROG_MAN + "/ionFinder/usage.txt";
	std::string const PROG_HELP_FILE = base::PROG_MAN + "/ionFinder/helpFile.roff";
	std::string const DEFAULT_FILTER_FILE_NAME = "DTASelect-filter.txt";
	
	//!default residues which are isobaric with a modification
	std::string const DEFAULT_AMBIGIOUS_RESIDUES = "";
	std::string const CIT_AMB_RESIDUES = "NQ";
	
	double const CIT_NL_MASS = 43.0058;
	double const DEFAULT_NEUTRAL_LOSS_MASS = CIT_NL_MASS;
	
	class Params;
	
	class Params : public base::ParamsBase{
	public:
		typedef std::map<std::string, std::string> FilterFilesType;
		
	private:
		std::string _parentDir;
		FilterFilesType _filterFiles;
		//!Default name of DTAFilter filter file to search for
		std::string _dtaFilterBase;
		//! mass of neutral loss to search for
		double _neutralLossMass;
		std::string _ambigiousResidues;
		
		//! should reverse peptide matches be considered
		bool _includeReverse;
		//!Should annotaed spectra be printed?
		bool _printSpectraFiles;
		//!Should NL ions be search for?
		bool _calcNL;
		
		//! was a input directory specified
		bool _inDirSpecified;
		
		//! names of folders to read
		std::vector<std::string> _inDirs;
		
		//! number of thread to use. Default is std::thread::hardware_concurrency() / 2
		unsigned int _numThread;
		
		//! path of fasta file to get modified residue numbers
		std::string _fastaFile;
		
		bool getFlist();
		unsigned int computeThreads() const;
		
	public:
		
		Params() : ParamsBase(PROG_USAGE_FNAME, PROG_HELP_FILE){
			_parentDir = "";
			_fastaFile = "";
			_includeReverse = false;
			_printSpectraFiles = false;
			_calcNL = true;
			_dtaFilterBase = DEFAULT_FILTER_FILE_NAME;
			_neutralLossMass = DEFAULT_NEUTRAL_LOSS_MASS;
			_ambigiousResidues = DEFAULT_AMBIGIOUS_RESIDUES;
			ofname = "peptide_cit_stats.tsv";
			_numThread = computeThreads();
			_inDirSpecified = false;
		}
		
		//modifers
		bool getArgs(int, const char* const[]);
		
		//properties
		const FilterFilesType& getFilterFiles() const{
			return _filterFiles;
		}
		bool getIncludeReverse() const{
			return _includeReverse;
		}
		bool getCalcNL() const{
			return _calcNL;
		}
		std::string getInputModeIndependentParentDir() const;
		double getNeutralLossMass() const{
			return _neutralLossMass;
		}
		std::string getAmbigiousResidues() const{
			return _ambigiousResidues;
		}
		std::string makeOfname() const{
			if(_inDirSpecified)
				return _wd + "/" + ofname;
			else{
				assert(_inDirs.size() == 1);
				return _inDirs.back() + "/" + ofname;
			}
		}
		bool getPrintSpectraFiles() const{
			return _printSpectraFiles;
		}
		unsigned int getNumThreads() const{
			return _numThread;
		}
		bool getInDirSpecified() const{
			return _inDirSpecified;
		}
		std::string getFastaFile() const{
			return _fastaFile;
		}
	};
}

#endif /* params_hpp */
