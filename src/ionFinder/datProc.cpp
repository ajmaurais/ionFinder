//
//  datProc.cpp
//  citFinder
//
//  Created by Aaron Maurais on 12/10/18.
//  Copyright © 2018 Aaron Maurais. All rights reserved.
//

#include <ionFinder/datProc.hpp>

/**
 Clear PeptideFragmentsMap.
 */
void IonFinder::PeptideFragmentsMap::clear()
{
	fragmentMap.clear();
	_sequence.clear();
}

/**
 \brief Generate fragment map from sequence.
 
 From a full sequence given, the sequence of all b and y ions are calculated
 and stored in the object.
 \param sequence the peptide sequence to calculate fragments from
 */
void IonFinder::PeptideFragmentsMap::populateMap(std::string sequence)
{
	_sequence = sequence;
	size_t len = _sequence.length();
	for(int i = 0; i < len; i++)
	{
		std::string beg = _sequence.substr(0, i + 1);
		std::string end = _sequence.substr(i);
		
		fragmentMap["b" + std::to_string(i+1)] = beg;
		
		if(i == 0)fragmentMap["M"] = end; //if first y ion, y ion is actually M
		else fragmentMap["y" + std::to_string(len-i)] = end; //regular y ion
	}
}

std::string IonFinder::PeptideFragmentsMap::getIonSeq(std::string searchStr) const{
	return fragmentMap.at(searchStr);
}

std::string IonFinder::PeptideFragmentsMap::getIonSeq(char b_y, int num) const{
	std::string searchStr = std::string(1, b_y) +
		(b_y == 'M' || b_y == 'm' ? "" : std::to_string(num)); //only add num if not M ion
	return getIonSeq(searchStr);
}

void IonFinder::RichFragmentIon::calcSequence(const PeptideFragmentsMap& pepMap){
	sequence = pepMap.getIonSeq(b_y, num);
}

/**
 \brief Align ref sequence to query sequence.
 
 Two peptide sequences as strings are aligned.
 \param ref The reference sequence to search.
 \param query The sequence to search for.
 \param beg Beginning of match. If no match, left unchanged.
 \param end End of match. If no match, left unchanged.
 \return bool representing if match was found.
 */
bool IonFinder::allignSeq(const std::string& ref, const std::string& query, size_t& beg, size_t& end)
{
	//find query in match
	size_t match = ref.find(query);
	if(match == std::string::npos) return false;
	
	beg = match;
	end = match + query.length() - 1;
	
	return true;
}

void IonFinder::PeptideStats::initStats()
{
	for(IonType i = IonType::First; i != IonType::Last; ++i){
		ionTypesCount[i] = IonTypeDatType("", 0);
	}
	
	containsCit = "false";
}

void IonFinder::PeptideStats::initModLocs(const char* diffmods)
{
	size_t modLoc = std::string::npos;
	bool modFound = false;
	size_t len = fullSequence.length();
	std::string tempSeq = fullSequence;
	for(size_t i = 0; i < len; i++)
	{
		//check that current char is not a mod char
		for(const char* p = diffmods; *p; p++)
			if(tempSeq[i] == *p)
				throw std::runtime_error("Invalid peptide sequence!");
		
		//if not at second to last AA, search for diff mod
		if((i + 1) < tempSeq.length())
		{
			//iterate through diffmods
			for(const char* p = diffmods; *p; p++)
			{
				//check if next char in sequence is diff mod char
				if(tempSeq[i + 1] == *p)
				{
					modFound = true;
					modLoc = i;
					tempSeq.erase(i + 1, 1);
					break;
				}//end of if
			}//end of for
		}//end of if
		if(modFound) modLocs.push_back(int(modLoc));
		modFound = false;
	}//end of for
}

void IonFinder::PeptideStats::addChar(std::string toAdd, std::string& s, std::string fragDelim)
{
	if(s.empty())
		s = toAdd;
	else s += fragDelim + toAdd;
}

/**
 Add ionStr and increment ion count for ion
 @param ionStr ion string to add
 @param ion PeptideStats::IonTypeDatType to increment
 @param inc amt to add to ion count
 */
void IonFinder::PeptideStats::incrementIonCount(std::string ionStr,
												PeptideStats::IonTypeDatType& ion,
												int inc){
	addChar(ionStr, ion.first, _fragDelim);
	ion.second += inc;
}

/**
 Add modified residue to PeptideStats::modResidues
 \param mod Modified residue in the form <residue><number> as in C57
 */
void IonFinder::PeptideStats::addMod(std::string mod){
	addChar(mod, modResidues, _fragDelim);
}

/**
 Tests whether sequence contains ambiguous residues.
 @param ambResidues ambiguous residues to search for
 @param fragSeq sequence of fragment to search for ambiguous residues.
 @return True if an ambiguous residue is found.
*/
bool IonFinder::PeptideStats::containsAmbResidues(const std::string& ambResidues,
												  std::string fragSeq) const
{
	size_t len = fragSeq.length();
	for(int i = 0; i < len; i++)
		for(int j = 0; j < ambResidues.length(); j++)
			if(fragSeq[i] == ambResidues[j])
				return true;
	return false;
}

/**
 Add fragment sequence to PeptideStats.
 @param seq fragment ion to add
 @param ambResidues ambiguous residues to search for.
 */
void IonFinder::PeptideStats::addSeq(const IonFinder::RichFragmentIon& seq,
									 const std::string& ambResidues)
{
	//first check if seq is found in *this sequence
	size_t beg, end;
	if(!IonFinder::allignSeq(sequence, seq.getSequence(), beg, end))
		return;
	
	//increment total fragment ions found
	std::string ionStr = seq.getLabel(true);
	incrementIonCount(ionStr, ionTypesCount[IonType::FRAG]);
	
	//check if seq spans any modified residues
	for(auto it = modLocs.begin(); it != modLocs.end(); ++it)
	{
		//check if in span
		if(utils::inSpan(beg, end, *it))
		{
			//check if NL
			if(seq.isNL()){
				incrementIonCount(ionStr, ionTypesCount[IonType::DET_NL_FRAG]);
			}
			else
			{
				if(containsAmbResidues(ambResidues, seq.getSequence())){ //is ambModFrag
					incrementIonCount(ionStr, ionTypesCount[IonType::AMB_MOD_FRAG]);
				}
				else{ //is detFrag
					incrementIonCount(ionStr, ionTypesCount[IonType::DET_FRAG]);
				}
			}
		}
		else{
			if(seq.isNL()){ //is artifact NL frag
				incrementIonCount(ionStr, ionTypesCount[IonType::ART_NL_FRAG]);
			}
			else{ //is amg frag
				incrementIonCount(ionStr, ionTypesCount[IonType::AMB_FRAG]);
			}//end of else
		}//end of else
	}//end of for
}//end of fxn

/**
 
 @param peptides vector of peptides to analyze
 */
bool IonFinder::analyzeSequences(std::vector<Dtafilter::Scan>& scans,
								 const std::vector<PeptideNamespace::Peptide>& peptides,
								 std::vector<PeptideStats>& peptideStats,
								 const IonFinder::Params& pars)
{
	IonFinder::PeptideFragmentsMap fragmentMap;
	bool allSucess = true;
	bool addModResidues = pars.getFastaFile() != "";
	fastaFile::FastaFile seqFile;
	int nSeqNotFound = 0;
	if(addModResidues){
		std::cout << "\nReading FASTA file...";
		if(!seqFile.read(pars.getFastaFile())) return false;
		std::cout << "Done!" << NEW_LINE;
	}
	
	for(auto it = peptides.begin(); it != peptides.end(); ++it)
	{
		/*if(it->getSequence() == "RVMGPDFG"){
			std::cout << "Found!" << NEW_LINE;
			it->printFragments(std::cout);
		}*/
		
		fragmentMap.clear();
		fragmentMap.populateMap(it->getSequence());
		//it->printFragments(std::cout);
		
		//initilize new pepStat object
		IonFinder::PeptideStats pepStat(*it); //init pepStat
		pepStat._scan = &scans[it - peptides.begin()]; //add pointer to scan
		size_t nFragments = it->getNumFragments();
		
		//iterate throught ion fragmetns
		for(size_t i = 0; i < nFragments; i++)
		{
			/*if(it->getFragment(i).getIonStr(false) == "b14-43")
				std::cout << "Found!" << NEW_LINE;*/
			
			//skip if not found
			if(it->getFragment(i).getFound())
			{
				IonFinder::RichFragmentIon fragTemp(it->getFragment(i));
				fragTemp.calcSequence(fragmentMap);
				pepStat.addSeq(fragTemp, pars.getAmbigiousResidues());
			} //end of if
		}//end of for i
		pepStat.calcContainsCit();
		
		//add modified residues if fasta file was specified
		if(addModResidues)
		{
			for(auto it = pepStat.modLocs.begin(); it != pepStat.modLocs.end(); ++it)
			{
				bool found; //set to true if peptide and prot sequences are found in FastaFile
				std::string modTemp = seqFile.getModifiedResidue(pepStat._scan->getParentID(),
																 pepStat.sequence, int(*it), found);
				pepStat.addMod(modTemp);
				if(!found)
				{
					nSeqNotFound++;
					if(pars.getVerbose())
						std::cerr << NEW_LINE;
				}
			}
		}
		peptideStats.push_back(pepStat);
	}//end if for it
	if(nSeqNotFound > 0){
		std::cerr << NEW_LINE << nSeqNotFound << " protein sequences not found in " <<
		pars.getFastaFile() << NEW_LINE;
	}
		
	
	return allSucess;
}//end of fxn

/**
 Search parent ms2 files in \p scans for predicted fragment ions. <br><br>
 Analysis is performed in parallel in number of threads in Params::_numThread. <br>
 \p scans is split up evenly across each thread.
 
 \param scans populated list of identified ms2 scans to search for
 \param peptides empty list of peptides to annotate
 \param pars Params object for information on how to perform analysis
 \return true is all file I/O was sucessful.
 */
bool IonFinder::findFragmentsParallel(std::vector<Dtafilter::Scan>& scans,
									  std::vector<PeptideNamespace::Peptide>& peptides,
									  const IonFinder::Params& pars)
{
	unsigned int const nThread = pars.getNumThreads();
	size_t const nScans = scans.size();
	size_t peptidePerThread = nScans / nThread;
	if(nScans % nThread != 0)
		peptidePerThread += 1;
	std::atomic<size_t> scansIndex(0); //used to update progress for findFragmentsProgress
	
	//init threads
	std::vector<std::thread> threads;
	bool* sucsses = new bool[nThread];
	
	//read ms2s
	std::cout << "Reading parent ms2 files...";
	Ms2Map ms2Map;
	if(!IonFinder::readMs2s(ms2Map, scans)) return false;
	std::cout << "Done!\n";
	
	//split up input data for each thread
	std::vector<PeptideNamespace::Peptide>* splitPeptides = new std::vector<PeptideNamespace::Peptide>[nThread];
	size_t begNum, endNum ;
	unsigned int threadIndex = 0;
	for(size_t i = 0; i < nScans; i += peptidePerThread)
	{
		begNum = i;
		endNum = (begNum + peptidePerThread > nScans ? nScans : begNum + peptidePerThread);

		//spawn thread
		assert(threadIndex < nThread);
		splitPeptides[threadIndex] = std::vector<PeptideNamespace::Peptide>();
		threads.push_back(std::thread(IonFinder::findFragments_threadSafe, std::ref(scans), begNum, endNum,
									  ms2Map,
									  std::ref(splitPeptides[threadIndex]), std::ref(pars),
									  sucsses + threadIndex, std::ref(scansIndex)));
		threadIndex++;
	}
	
	//spawn progress function
	if(!pars.getVerbose())
		threads.push_back(std::thread(IonFinder::findFragmentsProgress, std::ref(scansIndex), nScans,
									  nThread, PROGRESS_SLEEP_TIME));
	
	//join threads
	for(auto it = threads.begin(); it != threads.end(); ++it){
		it->join();
	 }
	
	//concat split peptides into one vector
	peptides.clear();
	for(unsigned int i = 0; i < nThread; i++){
		if(!sucsses[i])
			return false;
		peptides.insert(peptides.end(), splitPeptides[i].begin(), splitPeptides[i].end());
	}
	
	delete [] splitPeptides;
	delete [] sucsses;
	return true;
}

/**
 Prints progress bar during findFragmentsParallel. <br>
 \p scansIndex is updated concurently by each thread.
 \param scansIndex number of scans completed
 \param count total number of scans to search for
 \param sleepTime time before next update is printed (in seconds)
 */
void IonFinder::findFragmentsProgress(std::atomic<size_t>& scansIndex, size_t count,
									  unsigned int nThread, int sleepTime)
{
	size_t lastIndex = scansIndex.load();
	size_t curIndex = lastIndex;
	int noChangeIterations = 0;
	
	std::cout << "\nSearching ms2s for fragment ions using " << nThread << " thread(s)...\n";
	while(scansIndex < count)
	{
		curIndex = scansIndex.load();
		
		if(lastIndex == curIndex)
			noChangeIterations++;
		else noChangeIterations = 0;
		
		if(noChangeIterations > IonFinder::MAX_PROGRESS_ITTERATIONS)
			return;
		
		utils::printProgress(float(curIndex) / float(count));
		std::this_thread::sleep_for(std::chrono::seconds(sleepTime));
		
		lastIndex = curIndex;
	}
	utils::printProgress(float(scansIndex.load()) / float(count));
	std::cout << NEW_LINE;
	std::cout << "Done!" << NEW_LINE;
}

/**
 Find peptide fragment ions in ms2 files.
 @param scans Populated vector of scan objects to search for
 @param peptides empty vector of peptides to be filled from data in scans.
 @param pars IonFinder params object.
 @return true if successful
 */
bool IonFinder::findFragments(std::vector<Dtafilter::Scan>& scans,
							  std::vector<PeptideNamespace::Peptide>& peptides,
							  IonFinder::Params& pars)
{
	bool* success = new bool(false);
	std::atomic<size_t> scansIndex;
	IonFinder::Ms2Map ms2Map;
	if(!readMs2s(ms2Map, scans)) return false;
	IonFinder::findFragments_threadSafe(scans, 0, scans.size(), ms2Map, peptides, pars,
										success, scansIndex);
	
	bool ret = *success;
	delete success;
	return ret;
}

bool IonFinder::readMs2s(IonFinder::Ms2Map& ms2Map,
						 const std::vector<Dtafilter::Scan>& scans)
{
	std::string curWD;
	
	//first get unique names of ms2 files to read
	size_t len = scans.size();
	std::vector<std::string> fileNamesList(len);
	for(size_t i = 0; i < len; i++)
		fileNamesList[i] = scans[i].getPrecursorFile();
	auto it = std::unique(fileNamesList.begin(), fileNamesList.end());
	fileNamesList.resize(std::distance(fileNamesList.begin(), it));
	
	//read ms2 files
	ms2Map.clear();
	for(auto it = fileNamesList.begin(); it != fileNamesList.end(); ++it)
	{
		ms2Map[*it] = ms2::Ms2File();
		if(!ms2Map[*it].read(*it)){
			std::cerr << "Error reading ms2 files!" << NEW_LINE;
			return false;
		}
	}
	return true;
}

/**
 Find peptide fragment ions in ms2 files.
 Function is thread safe.
 @param peptides empty vector of peptides to be filled from data in scans.
 @param beg index of beginning of scan vector
 @param end index of end of scan vector
 @param pars IonFinder params object.
 @param success set to true if function was successful
 */
void IonFinder::findFragments_threadSafe(std::vector<Dtafilter::Scan>& scans,
										 const size_t beg, const size_t end,
										 const IonFinder::Ms2Map& ms2Map,
										 std::vector<PeptideNamespace::Peptide>& peptides,
										 const IonFinder::Params& pars,
										 bool* success, std::atomic<size_t>& scansIndex)
{
	*success = false;
	std::string curSample;
	std::string curWD;
	std::string spFname;
	aaDB::AADB aminoAcidMasses;
	ms2::Spectrum spectrum;
	
	for(size_t i = beg; i < end; i++)
	{
		//first get current wd name
		curWD = utils::dirName(scans[i].getPrecursorFile());
		spFname = curWD + "/sequest.params";
		
		if(curSample != scans[i].getSampleName())
		{
			//re-init Peptide::AminoAcidMasses for each sample
			curWD = utils::dirName(scans[i].getPrecursorFile());
			spFname = curWD + "/sequest.params";
			
			//read sequest params file and init aadb
			PeptideNamespace::initAminoAcidsMasses(pars, spFname, aminoAcidMasses);
		}//end if
		curSample = scans[i].getSampleName();
		
		//initialize peptide object for current scan
		peptides.push_back(PeptideNamespace::Peptide(scans[i].getSequence()));
		peptides.back().initialize(pars, aminoAcidMasses);
		
		if(pars.getCalcNL()){
			//calculate neutral loss combinations
			int nMods = peptides.back().getNumMod();
			double nlMass = pars.getNeutralLossMass();
			std::vector<double> neutralLossIons;
			for(int i = 1; i <= nMods; i++)
				neutralLossIons.push_back(i * nlMass);
			
			//add neutral loss fragments to current peptide
			peptides.back().addNeutralLoss(neutralLossIons);
		}
		
		//load spectrum
		auto ms2FileIt = ms2Map.find(scans[i].getPrecursorFile());
		if(ms2FileIt == ms2Map.end()){
			throw std::runtime_error("Key error in Ms2Map!");
			return;
		}
		if(!ms2FileIt->second.getScan(scans[i].getScanNum(), spectrum)){
			throw std::runtime_error("Error reading scan!");
			return;
		}
		scans[i].setPrecursorMZ(spectrum.getPrecursorMZ());
		scans[i].setPrecursorScan(spectrum.getPrecursorScan());
		
		//spectrum.labelSpectrum(peptides.back(), pars, true); //removes unlabeled ions from peptide
		spectrum.labelSpectrum(peptides.back(), pars);
		
		//print spectra file
		if(pars.getPrintSpectraFiles())
		{
			std::string dirNameTemp = (pars.getInDirSpecified() ? pars.getWD() : curWD) + "/spectraFiles";
			if(!utils::dirExists(dirNameTemp))
				if(!utils::mkdir(dirNameTemp.c_str())){
					throw std::runtime_error("Failed to make dir: " + dirNameTemp);
					return;
				}
			spectrum.normalizeIonInts(100);
			spectrum.calcLabelPos();
			spectrum.setCharge(scans[i].getCharge());
			
			std::string temp = dirNameTemp + "/" + utils::baseName(scans[i].getOfname());
			std::ofstream outF((temp).c_str());
			if(!outF){
				throw std::runtime_error("Failed to write spectrum!");
				return;
			}
			spectrum.printLabeledSpectrum(outF, true);
		}
		scansIndex++;
	} //end of for
	
	*success = true;
	return;
}

void IonFinder::PeptideStats::calcContainsCit()
{
	containsCit = "false";
	
	//is the peptide modified?
	if(modLocs.size() == 0) return;
	
	//is c terminal most cit modification on the c terminus?
	//modLocs.back() works because modLocs are added in the order
	//they appear in the sequence
	if(modLocs.back() == sequence.length() - 1) return;
	
	//is there more than 1 determining NLs?
	if(ionTypesCount[IonType::DET_NL_FRAG].second > 1){
		containsCit = "true";
		return;
	}
	
	//are there 1 or more determining NLs or determining frags?
	if(ionTypesCount[IonType::DET_NL_FRAG].second >= 1 ||
	   ionTypesCount[IonType::DET_FRAG].second >= 1){
		containsCit = "likely";
		return;
	}
	
	//are there 1 more ambiguous fragments?
	if(ionTypesCount[IonType::AMB_FRAG].second >= 1){
		containsCit = "ambiguous";
		return;
	}
}

std::string IonFinder::PeptideStats::ionTypeToStr(const IonType& it)
{
	switch(it){
		case IonType::FRAG: return ION_TYPES_STR[0];
			break;
		case IonType::DET_FRAG: return ION_TYPES_STR[1];
			break;
		case IonType::AMB_MOD_FRAG: return ION_TYPES_STR[2];
			break;
		case IonType::DET_NL_FRAG: return ION_TYPES_STR[3];
			break;
		case IonType::AMB_FRAG: return ION_TYPES_STR[4];
			break;
		case IonType::ART_NL_FRAG: return ION_TYPES_STR[5];
			break;
		case IonType::Last: return "Last";
			break;
	}
}

/**
 Prints peptide stats to out.
 @param stats Peptide stats to print.
 @param pars initalized IonFinder::Params object
 */
bool IonFinder::printPeptideStats(const std::vector<PeptideStats>& stats,
								  const IonFinder::Params& pars)
{
	//assert(out);
	std::ofstream outF (pars.makeOfname());
	if(!outF) return false;
	
	typedef IonFinder::PeptideStats::IonType itcType;
	//build stat names vector
	std::vector<std::string> statNames;
	
	if(pars.getCalcNL())
		statNames.push_back("contains_Cit");
	else statNames.push_back("contains_mod");
	
	//determine when to stop printing peptide stats based on analysis performed
	std::vector<itcType> _pepStats; //used to store relavent peptide stats based on params
	//defaults
	_pepStats.push_back(itcType::FRAG);
	_pepStats.push_back(itcType::DET_FRAG);
	_pepStats.push_back(itcType::AMB_FRAG);
	//conditional stats
	if(!pars.getAmbigiousResidues().empty())
		_pepStats.push_back(itcType::AMB_MOD_FRAG);
	if(pars.getCalcNL()){
		_pepStats.push_back(itcType::DET_NL_FRAG);
		_pepStats.push_back(itcType::ART_NL_FRAG);
	}
	
	//append peptide stats names to headers
	int statLen = 0;
	for(auto it = _pepStats.begin(); it != _pepStats.end(); ++it){
		statNames.push_back("n" +
							std::string(1, (char)std::toupper(ION_TYPES_STR[utils::as_integer(*it)][0])) +
							ION_TYPES_STR[utils::as_integer(*it)].substr(1));
		statLen++;
	}
	statNames.insert(statNames.end(), ION_TYPES_STR, ION_TYPES_STR + statLen);
	
	std::string otherHeaders = "protein_ID parent_protein protein_description full_sequence sequence parent_mz is_modified modified_residue charge unique xCorr scan precursor_scan parent_file sample_name";
	std::vector<std::string> oHeaders;
	utils::split(otherHeaders, ' ', oHeaders);
	std::vector<std::string> headers;
	headers.insert(headers.end(), oHeaders.begin(), oHeaders.end());
	headers.insert(headers.end(), statNames.begin(), statNames.end());
	
	//print headers
	size_t len = headers.size();
	for(size_t i = 0; i < len; i++){
		if(i == 0)
			outF << headers[i];
		else outF << OUT_DELIM << headers[i];
	}
	outF << NEW_LINE;
	
	//print data
	for(auto it = stats.begin(); it != stats.end(); ++it)
	{
		//scan data
		outF << it->_scan->getParentID() <<
		OUT_DELIM << it->_scan->getParentProtein() <<
		OUT_DELIM << it->_scan->getParentDescription() <<
		OUT_DELIM << it->_scan->getFullSequence() <<
		OUT_DELIM << it->_scan->getSequence() <<
		OUT_DELIM << it->_scan->getPrecursorMZ() <<
		OUT_DELIM << (it->modLocs.size() > 0) <<
		OUT_DELIM << it->modResidues <<
		OUT_DELIM << it->_scan->getCharge() <<
		OUT_DELIM << it->_scan->getUnique() <<
		OUT_DELIM << it->_scan->getXcorr() <<
		OUT_DELIM << it->_scan->getScanNum() <<
		OUT_DELIM << it->_scan->getPrecursorScan() <<
		OUT_DELIM << utils::baseName(it->_scan->getPrecursorFile()) <<
		OUT_DELIM << it->_scan->getSampleName();
		
		//peptide analysis data
		outF << OUT_DELIM;
		if(pars.getCalcNL())
			 outF << it->containsCit;
		else{
			outF << (it->ionTypesCount.at(itcType::DET_FRAG).second > 0);
		}
		
		for(auto it2 = _pepStats.begin(); it2 != _pepStats.end(); ++it2)
			outF << OUT_DELIM << it->ionTypesCount.at(*it2).second;
		
		for(auto it2 = _pepStats.begin(); it2 != _pepStats.end(); ++it2)
			outF << OUT_DELIM << it->ionTypesCount.at(*it2).first;
		
		outF << NEW_LINE;
	}
	return true;
}
