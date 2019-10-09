//
//  peptide.hpp
//  ms2_anotator
//
//  Created by Aaron Maurais on 3/20/17.
//  Copyright © 2017 Aaron Maurais. All rights reserved.
//

#ifndef peptide_hpp
#define peptide_hpp

#include <map>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <string>
#include <iomanip>

#include <utils.hpp>
#include <aaDB.hpp>
#include <paramsBase.hpp>
#include <sequestParams.hpp>

namespace PeptideNamespace{
	
	//forward class declarations
	class Species;
	class Ion;
	class AminoAcid;
	class Peptide;
	class FragmentIon;
	
	const double H_MASS = 1.00732;
	//const double H_MASS = 1.00783;
	
	typedef std::vector<AminoAcid> PepIonVecType;
	typedef PepIonVecType::const_iterator PepIonIt;
	
	//forward function declarations
	//peptide mass functions
	double calcMass(double mz, int charge);
	double calcMZ(double mass, int charge);
	double calcMass(std::string sequence);
	double calcMass(PepIonIt begin, PepIonIt end);
	
	//Modification helper functions
	std::string concatMods(PepIonIt begin, PepIonIt end);
	
	//mass calculation functions
	void initAminoAcidsMasses(const base::ParamsBase& pars, std::string seqParFname, aaDB::AADB&);
	void initAminoAcidsMasses(const base::ParamsBase&, aaDB::AADB&);
	
	//class declarations
	
	//!base class for peptide species
	class Species{
	protected:
		double mass;
		
	public:
		Species(){
			mass = 0;
		}
		
		//properties
		double getMass() const{
			return mass;
		}
	};
	
	//!base class for all ions
	class Ion : public Species{
	protected:
		int charge;
	public:
		//constructors
		Ion() : Species(){
			charge = 0;
		}
		~Ion(){}
		
		//modifers
		void initializeFromMZ(double _mz, int _charge){
			mass = calcMass(_mz, _charge);
			charge = _charge;
		}
		void initalizeFromMass(double _mass, int _charge = 1){
			mass = _mass;
			charge = _charge;
		}
		
		double getMZ(int _charge) const{
			return calcMZ(mass, _charge);
		}
		double getMZ() const{
			return getMZ(charge);
		}
		int getCharge() const{
			return charge;
		}
		std::string makeChargeLable() const;
	};//end of class
	
	class AminoAcid : public Ion{
	protected:
		//!Does the amino acid bear a static modification?
		bool _staticMod;
		//!Does amino acid have a dynamic modification?
		bool _dynamicMod;
		//!str used to represent mod
		char _mod;
		//!modification mass change (can be potitive or negative)
		double _modMass;
		
		void _addMod(double modMass);
	public:
		AminoAcid() : Ion() {
			_staticMod = false;
			_dynamicMod = false;
			_modMass = 0;
			_mod = '\0';
		}
		AminoAcid(double mass) : Ion() {
			initalizeFromMass(mass);
			_staticMod = false;
			_dynamicMod = false;
			_modMass = 0;
			_mod = '\0';
		}
		
		void setDynamicMod(char mod, double modMass);
		void addStaticMod(double modMass);
		
		//std::string makeModLable() const;
		double getModMass() const{
			return _modMass;
		}
		//!Get base mass of amino acid + modMass
		double getTotalMass() const{
			return _modMass + mass;
		}
		char getMod() const{
			return _mod;
		}
		//!Does AA have any modification (static or dynamic)
		bool isModified() const{
			return _staticMod || _dynamicMod;
		}
		//!Does AA have dynamic modification
		bool hasDynamicMod() const{
			return _dynamicMod;
		}
		
	};
	
	//!Used to represent b and y peptide ions
	class FragmentIon : public Ion{
	public:
		enum class IonType{B, Y, M, B_NL, Y_NL, M_NL};
		
	protected:
		char _b_y;
		int _num;
		//!Represents all symbol for modification 
		 /** i.e. if two modifications are present, mod will be "**" */
		std::string _mod;
		//!Was peptide fragment found in ms2 spectra?
		bool _found;
		IonType _ionType;
		 //!Represents neutral loss mass
		double _nlMass;
		//!Represents multiples of base neutral loss mass on peptide
		size_t _numNl;
		//!Should ion label be included in spectrum?
		bool _includeLabel;
		
		//!sequence of fragment
		std::string _sequence;
		//!index of beginning of fragment relative to full sequence
		size_t _beg;
		//!index of end of fragment relative to full sequence
		size_t _end;
		
		void _initFragSpan(const std::string&);
		
	public:
		//!blank constructor
		FragmentIon() : Ion(){
			_b_y = '\0';
			_num = 0;
			_mod = "";
			_found = false;
			_ionType = IonType::B;
			_nlMass = 0.0;
			_numNl = 0;
			_sequence = "";
			_beg = std::string::npos;
			_end = std::string::npos;
			_includeLabel = false;
		}
		FragmentIon(char b_y, int num, int charge, double mass,
					std::string mod, const std::string& pepSequence) : Ion() {
			_b_y = b_y;
			_num = num;
			_mod = mod;
			_nlMass = 0;
			_numNl = 0;
			initalizeFromMass(mass, charge);
			_found = false;
			_ionType = strToIonType(b_y);
			_initFragSpan(pepSequence);
			_includeLabel = true;
		}
		//!copy constructor
		FragmentIon(const FragmentIon& rhs){
			_b_y = rhs._b_y;
			_num = rhs._num;
			_mod = rhs._mod;
			_found = rhs._found;
			_ionType = rhs._ionType;
			_nlMass = rhs._nlMass;
			_numNl = rhs._numNl;
			charge = rhs.charge;
			mass = rhs.mass;
			_sequence = rhs._sequence;
			_beg = rhs._beg;
			_end = rhs._end;
			_includeLabel = rhs._includeLabel;
		}
		~FragmentIon() {}
		
		void setFound(bool boo){
			_found = boo;
		}
		void setIonType(IonType it){
			_ionType = it;
		}
		void setForceLabel(bool boo){
			_includeLabel = boo;
		}
		
		//properties
		double getMZ() const{
			if(_b_y == 'b')
				return (mass + ((charge - 1) * PeptideNamespace::H_MASS)) / charge;
			return Ion::getMZ(charge);
		}
		std::string getLabel(bool includeMod = true, std::string chargeSep = " ") const;
		std::string getFormatedLabel() const;
		char getBY() const{
			return _b_y;
		}
		//!Get fragment ion number
		int getNum() const{
			return _num;
		}
		//!Get copy of FragmentIon::mod
		std::string getMod() const{
			return _mod;
		}
		//!Get number of modifications on fragment
		size_t getNumMod() const{
			return _mod.length();
		}
		//!Get number of neutral loss multiples on fragment
		size_t getNumNl() const{
			return _numNl;
		}
		bool getFound() const{
			return _found;
		}
		IonType getIonType() const{
			return _ionType;
		}
		bool getIncludeLabel() const{
			return _includeLabel;
		}
		static IonType strToIonType(std::string);
		static IonType strToIonType(char c){
			return strToIonType(std::string(1, c));
		}
		std::string ionTypeToStr() const;
		std::string getNLStr() const;
		bool isModified() const{
			return _mod == "";
		}
		//!returns true if fragment is neutral loss ion
		bool isNL() const{
			return _ionType == IonType::B_NL ||
			_ionType == IonType::Y_NL ||
			_ionType == IonType::M_NL;
		}
		//!\return true if fragment is parent ion or parent neutral loss
		bool isM() const{
			return _ionType == IonType::M || _ionType == IonType::M_NL;
		}
		FragmentIon makeNLFrag(double lossMass, size_t numNL) const;
		
		std::string getSequence() const{
			return _sequence;
		}
		size_t getBegin() const{
			return _beg;
		}
		size_t getEnd() const{
			return _end;
		}
	};
	
	//!Used to store fragment data for each peptide.
	class Peptide : public Ion{
	private:
		typedef std::vector<FragmentIon> FragmentIonType;
		typedef FragmentIonType::const_iterator FragmentIonItType;
		
		std::string sequence;
		std::string fullSequence;
		std::vector<AminoAcid> aminoAcids;
		bool initialized;
		FragmentIonType fragments;
		//!number of modified residues
		int nMod;
		//!Locations of dynamic modifications on peptide sequence
		std::vector<size_t> modLocs;
		
		double parseStaticMod(size_t);
		void fixDiffMod(const aaDB::AADB& aminoAcidsMasses,
						const char* diffmods = "*");
		int nModsInSpan(size_t beg, size_t end) const;
	public:
		//constructors
		Peptide() : Ion(){
			sequence = "";
			fullSequence = sequence;
			initialized = false;
			nMod = 0;
		}
		Peptide(std::string _sequence) : Ion(){
			sequence = _sequence;
			fullSequence = sequence;
			initialized = false;
			nMod = 0;
		}
		~Peptide() {}
		
		//modifers
		void initialize(const base::ParamsBase&, const aaDB::AADB& aadb,
						bool _calcFragments = true);
		void calcFragments(int minCharge, int maxCharge,
						   const aaDB::AADB& aminoAcidsMasses);
		void addNeutralLoss(double losses, bool labelDecoyNL = false);
		double calcMass(const aaDB::AADB& aminoAcidsMasses);
		void printFragments(std::ostream&) const;
		
		void setFound(size_t i, bool boo){
			fragments[i].setFound(boo);
		}
		void removeUnlabeledFrags();
		
		//properties
		std::string getSequence() const{
			return sequence;
		}
		std::string getFullSequence() const{
			return fullSequence;
		}
		size_t getNumFragments() const{
			return fragments.size();
		}
		double getFragmentMZ(size_t i) const{
			return fragments[i].getMZ();
		}
		std::string getFragmentLabel(size_t i) const{
			return fragments[i].getLabel();
		}
		std::string getFormatedLabel(size_t i) const{
			return fragments[i].getFormatedLabel();
		}
		bool getIncludeLabel(size_t i) const{
			return fragments[i].getIncludeLabel();
		}
		template<typename _Tp>
		std::string getFormatedLabel(size_t i, _Tp num) const{
			return fragments[i].getFormatedLabel(num);
		}
		char getBY(size_t i) const{
			return fragments[i].getBY();
		}
		bool getFound(size_t i) const{
			return fragments[i].getFound();
		}
		const FragmentIon& getFragment(size_t i) const{
			return fragments[i];
		}
		int getNumMod() const{
			return nMod;
		}
		const std::vector<size_t>& getModLocs() const{
			return modLocs;
		}
	};//end of class
	
}//end of namespace

#endif /* peptide_hpp */
