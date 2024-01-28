#ifndef _CALM_H_
#define _CALM_H_

#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <TKey.h>
#include "TRandom2.h"
#include "TGenPhaseSpace.h"
#include "ParticleType.h"
#include "ParticleDB.h"
#include "Particle.h"
#include "THGlobal.h"
#include "TF1.h"
#include "TH1F.h"
#include "TFile.h"
//--------reggae
#include "reggae.h"
#include "specrel.h"
#include "ConfigurationHolder.h"

using namespace std;

enum eEventType
{
	GLOBAL,
	MINIJETS_GLOBAL,
	MINIJETS_LOCAL,
	GLOBAL_REGGAE,
	MINIJETS_GLOBAL_REGGAE,
	MINIJETS_LOCAL_REGGAE
};

class CALM
{
public:
	CALM();
	~CALM();
	int GenerateParticles(ParticleDB *aPartDB, int aMultBinMin, int aMultBinMax, double aEnergy, std::list<Particle> *aParticles, eEventType aEventType = GLOBAL, int aJetNumber = 2);
	//void SetMultiplicities(ParticleDB *aDB, int aMultBinMin, int aMultBinMax);
	//void Randomize();

private:
	TRandom2 *mRandom; ///< Random number generator
	// values constant (what can be generated etc)
	int mNpart;				  ///< Number of kinds of particles (pions, kaons, nuclides etc.)
	double *mNmean;			  ///< Mean values for the particle yields dN/dy of each particle count -> taken from ALICE
	double mRapidityInterval; ///< Interval of rapidity
	double *mXYZ;			  ///< To generate distance from primary vertex
	int *mNpartkinds;		  ///< Number of particles for each kind
	string **mNames;		  ///< Names of particles to be generated
	
	TFile *singleEnergyFile;	  ///< Energy distribution of single particle
	TFile *pionsFile;	  ///< Custom distribution function for pions
	TFile *kaonsFile;      ///< Custom distribution function for kaons
	TFile *nucleonsFile;   ///< Custom distribution function for nucleons
	TFile *lambdasFile;    ///< Custom distribution function for lambdas
	
	TF1 *singleEnergyDistr;	  ///< Energy distribution of single particle
	TF1 *pionsMultDistr;	  ///< Custom distribution function for pions
	TF1 *kaonsMultDistr;      ///< Custom distribution function for kaons
	TF1 *nucleonsMultDistr;   ///< Custom distribution function for nucleons
	TF1 *lambdasMultDistr;    ///< Custom distribution function for lambdas
	
	double Etot; 			  ///< Total momentum to be distributed among particles
	// values for this event
	vector<string> mParticlesThisEvent; ///< List of particle names for one event
	ConfigurationHolder* eventConfig; ///< ConfigurationHolder with parameters from file
	int *(CALM::*GetMultiplicitiesOfPartciles)(int, int, int &); ///< Pointer to the function generating amount of particles: CustomMult or AlicePoissonMult

private:
	//one of these methods have to be assaigned into GetMultiplicitiesOfPartciles pointer
	int *CustomMult(int aMultBinMin, int aMultBinMax, int &Nsum);
	int *AlicePoissonMult(int aMultBinMin, int aMultBinMax, int &Nsum);

	int **GetTypesForParticles(int *Nrand, ParticleDB *aPartDB);
	double **GetVertexXYZ(int Nsum);
	void SetTotalEnergy(int Nsum, double aEnergy);

	double *GetMasses(int Nsum, ParticleDB *aPartDB);


	void SeparateJets(int Nsum, vector<double> *masses, vector<string> *names, ParticleDB *aPartDB, int MinMulti, int JetNumber);
	bool SeparateJets_LOCAL(int Nsum, vector<double> *masses, vector<string> *names, ParticleDB *aPartDB, int MinMulti, int JetNumber);
	
	bool TrySetEventDecay(int Nsum, double *masses, TGenPhaseSpace &event, double &TotEnergy);
	bool TrySetEventDecay_MINIJETS(int Nsum, vector<double> *masses, TGenPhaseSpace &event0, TGenPhaseSpace &event1, TGenPhaseSpace &event2, TGenPhaseSpace &event3, double &TotEnergy, double *divideEn, int JetNumber);

	bool FilterUnlikelyEvents(TGenPhaseSpace &event, double &weight);
	bool FilterUnlikelyEvents_MINIJETS(TGenPhaseSpace &event0, TGenPhaseSpace &event1, TGenPhaseSpace &event2, TGenPhaseSpace &event3, double &weight0, double &weight1, double &weight2, double &weight3, int JetNumber);

	bool ReggaeNegativeEnergyCheck(int Nsum, double *masses, vector4 en, vector4 *avec);
	bool ReggaeNegativeEnergyCheck_MINIJETS(int Nsum, vector<double> *masses, vector4 en, vector4 *avec0, vector4 *avec1, vector4 *avec2, vector4 *avec3, int NumberOfJets);

	void SaveAllParticles_GLOBAL(int Nsum, double weight, double **XYZrand, TGenPhaseSpace event, ParticleDB *aPartDB, list<Particle> *aParticles);
	void SaveAllParticles_MINIJETS(vector<double> *masses, vector<string> *names, double weight0, double weight1, double weight2, double weight3, double TotEnergy, double *divideEn, double **XYZrand, TGenPhaseSpace event0, TGenPhaseSpace event1, TGenPhaseSpace event2, TGenPhaseSpace event3, ParticleDB *aPartDB, list<Particle> *aParticles, eEventType aEventType, int JetNumber);
	void SaveAllParticles_GLOBAL_REGGAE(int Nsum, vector4 *avec, double **XYZrand, ParticleDB *aPartDB, list<Particle> *aParticles);
	void SaveAllParticles_MINIJETS_REGGAE(vector<double> *masses, vector<string> *names, vector4 *avec0, vector4 *avec1, vector4 *avec2, vector4 *avec3, double TotEnergy, double **XYZrand, ParticleDB *aPartDB, list<Particle> *aParticles, int NumberOfJets);


};

#endif


/*! @file CALM.h
 * @brief Representation of the program, handles all aspects of simulation.
 */
/*! @class CALM
 * @brief Representation of the program, handles all aspects of simulation.
 * 
 *
 * @fn CALM::CALM()
 * @brief Default constructor.
 *
 *
 * @fn CALM::~CALM()
 * @brief Destructor.
 *
 * @fn int CALM::GenerateParticles(ParticleDB *aPartDB, int aMultBinMin, int aMultBinMax, double aEnergy, list<Particle> *aParticles, eEventType aEventType = GLOBAL)
 * @brief Generates event and saves outcome particles. It returnes 0 if successed and 99 if not.
 * @param [in] aPartDB pointer to ParticleDB
 * @param [in] aMultBinMin minimum multiplicity of event
 * @param [in] aMultBinMax maxiumum multiplicity of event
 * @param [in] aEnergy maximum energy of event
 * @param [out] aParticles pointer to list<Particle> which contains all particle data generated in one event
 * @param [in] aEventType enum that specifies which CALM option will be performed
 * @retval 0 if simulation is correctly completed
 * @retval 99 if generated data has no predisposition to finish the event
 * 
 * 
 * 
 * @fn int *CALM::CustomMult(int aMultBinMin, int aMultBinMax, int &Nsum)
 * @brief Randomize multiciplity, from custom distributions, of each kind (pion, etc.) for event and returnes pointer to array of them
 * @param [in] aMultBinMin minimum multiplicity of event
 * @param [in] aMultBinMax maxiumum multiplicity of event
 * @param [out] Nsum reference to the variable holding amount of particles within event
 * 
 * @fn int *CALM::AlicePoissonMult(int aMultBinMin, int aMultBinMax, int &Nsum)
 * @brief Randomize multiciplity, from Poisson distributions, of each kind (pion, etc.) for event and returnes pointer to array of them
 * @param [in] aMultBinMin minimum multiplicity of event
 * @param [in] aMultBinMax maxiumum multiplicity of event
 * @param [out] Nsum reference to the variable holding amount of particles within event
 * 
 * 
 * 
 * @fn void CALM::GetTypesForParticles(int *Nrand, ParticleDB *aPartDB)
 * @brief Generates type for each particle and checks of the charge, strangeness and baryon number. Returnes pointer to array of pointers (particles) to the type
 * @param [in] Nrand pointer to array returned by CALM::GetMultiplicitiesOfPartciles()
 * @param [in] aPartDB pointer to particle data base
 * 
 * 
 * 
 * 
 * @fn double CALM::GetVertexXYZ(int Nsum)
 * @brief Generates XYZ vertex coordinates for each particle (from Gaussian distribution). Returnes array of arrays for each dimension.
 * @param [in] Nsum amount of all particles
 * 
 * 
 * 
 * @fn double CALM::SetTotalEnergy(int Nsum, double aEnergy)
 * @brief Generates energy for each particle (from eventConfig->singleEnergyDistr distribuation) and checks if the sum is not bigger then maximum energy
 * @param [in] Nsum amount of all particles
 * @param [in] aEnergy maximum energy of event
 * 
 * 
 * 
 * @fn double *CALM::GetMasses(int Nsum, ParticleDB *aPartDB)
 * @brief Returnes pointer to array filled by masses (read from aPartDB) of particles generated earler in the event
 * @param [in] Nsum amount of all particles
 * @param [in] aPartDB pointer to particle data base
 * 
 * 
 * 
 * @fn void CALM::SeparateJets(int Nsum, vector<double> *masses, vector<string> *names, ParticleDB *aPartDB)
 * @brief Randomly separate particles into two jets
 * @param [in] Nsum amount of all particles
 * @param [out] masses pointer to the array of vectors which will be filled by masses of particles for each jet
 * @param [out] names pointer to the array of vectors which will be filled by names of particles for each jet
 * @param [in] aPartDB pointer to particle data base
 * 
 * 
 *
 * @fn bool CALM::SeparateJets_LOCAL(int Nsum, vector<double> *masses, vector<string> *names, ParticleDB *aPartDB)
 * @brief Randomly separate particles into two jets and checks ConservAtion Laws for each jet. Returnes true if the check was successful. In case of 
 * 100 negative control results, method returns false and than CALM repeats whole event from the beginning.
 * @param [in] Nsum amount of all particles
 * @param [out] masses pointer to the array of vectors which will be filled by masses of particles for each jet
 * @param [out] names pointer to the array of vectors which will be filled by names of particles for each jet
 * @param [in] aPartDB pointer to particle data base
 * 
 * 
 * 
 * @fn bool CALM::TrySetEventDecay(int Nsum, double *masses, TGenPhaseSpace &event, double &TotEnergy)
 * @brief Tries to set decay of energy for event. In case of negative control result, CALM repeats whole event from the beginning.
 * @param [in] Nsum amount of all particles
 * @param [in] masses pointer to the array which is filled by masses of particles
 * @param [out] event reference to TGenPhaseSpace
 * @param [out] TotEnergy reference to variable holding total energy
 *  
 * 
 *
 * @fn bool CALM::TrySetEventDecay_MINIJETS(int Nsum, vector<double> *masses, TGenPhaseSpace &event0, TGenPhaseSpace &event1, double &TotEnergy, double *divideEn)
 * @brief Tries to set decay of energy for minijets. In case of negative control result, CALM repeats whole event from the beginning.
 * @param [in] Nsum amount of all particles
 * @param [in] masses pointer to the vector holding two arrays (jets) which are filled by masses of particles
 * @param [out] event0 reference to TGenPhaseSpace of first jet
 * @param [out] event1 reference to TGenPhaseSpace of second jet
 * @param [out] TotEnergy reference to variable holding total energy
 * @param [in] divideEn pointer to the array of multipliers used to set particles energy and jets boost energy
 * 
 * 
 * 
 * @fn bool CALM::FilterUnlikelyEvents(TGenPhaseSpace &event, double &weight)
 * @brief Filters the most unlikely events. If event has too small probability to happen, CALM repeats whole event from the beginning. 
 * @param [out] event reference to TGenPhaseSpace
 * @param [out] weight reference to variable holding weight of event
 * 
 * 
 * 
 * @fn bool CALM::FilterUnlikelyEvents_MINIJETS(TGenPhaseSpace &event0, TGenPhaseSpace &event1, double &weight0, double &weight1)
 * @brief Filters the most unlikely minijet events. If event has too small probability to happen, CALM repeats whole event from the beginning. 
 * @param [out] event0 reference to TGenPhaseSpace of first jet
 * @param [out] event1 reference to TGenPhaseSpace of second jet
 * @param [out] weight0 reference to variable holding weight of first jet
 * @param [out] weight1 reference to variable holding weight of second jet
 * 
 * 
 * 
 * @fn bool CALM::ReggaeNegativeEnergyCheck(int Nsum, double *masses, vector4 en, vector4 *avec)
 * @brief Sets decay of energy and checks negative energy. If negative energy still occures CALM repeats whole event from the beginning. 
 * @param [in] Nsum amount of all particles
 * @param [in] masses pointer to the array which is filled by masses of particles
 * @param [in] en energy vector
 * @param [out] avec pointer to array of 4D vectors holding energy for each partible
 * 
 * 
 * 
 * @fn bool CALM::ReggaeNegativeEnergyCheck_MINIJETS(int Nsum, vector<double> *masses, vector4 en, vector4 *avec0, vector4 *avec1)
 * @brief Sets decay of energy for each jet and checks negative energy. If negative energy still occures CALM repeats whole event from the beginning. 
 * @param [in] Nsum amount of all particles
 * @param [in] masses pointer to the vector holding two arrays (jets) which are filled by masses of particles
 * @param [in] en energy vector
 * @param [out] avec0 pointer to array of 4D vectors holding energy for each partible inside first jet
 * @param [out] avec1 pointer to array of 4D vectors holding energy for each partible inside second jet
 * 
 * 
 * 
 * @fn void CALM::SaveAllParticles_GLOBAL(int Nsum, double weight, double **XYZrand, TGenPhaseSpace event, ParticleDB *aPartDB, list<Particle> *aParticles)
 * @brief Pushes all particles into the array which will be saved to the file (for Genbod GLOBAL)
 * @param [in] Nsum amount of all particles
 * @param [in] weight variable holding weight of event
 * @param [in] XYZrand array holding XYZ vertex coordinates got from CALM::GetXYZ
 * @param [in] event event object
 * @param [in] aPartDB pointer to particle data base
 * @param [out] aParticles pointer to the list of data for all particles generated in current event
 * 
 * 
 * 
 * @fn void CALM::SaveAllParticles_MINIJETS(vector<double> *masses, vector<string> *names, double weight0, double weight1, double TotEnergy, double *divideEn,
 *  double **XYZrand, TGenPhaseSpace event0, TGenPhaseSpace event1, ParticleDB *aPartDB, list<Particle> *aParticles, eEventType aEventType)
 * @brief Pushes all particles into the array which will be saved to the file (for Genbod MINIJETS)
 * @param [in] masses pointer to the array of vectors which will be filled by masses of particles for each jet
 * @param [in] names pointer to the array of vectors which will be filled by names of particles for each jet
 * @param [in] weight0 reference to variable holding weight of first jet
 * @param [in] weight1 reference to variable holding weight of second jet
 * @param [in] TotEnergy variable holding total energy
 * @param [in] divideEn pointer to the array of multipliers used to set particles energy and jets boost energy
 * @param [in] XYZrand array holding XYZ vertex coordinates got from CALM::GetXYZ
 * @param [in] event0 event object for first jet
 * @param [in] event1 event object for second jet
 * @param [in] aPartDB pointer to particle data base
 * @param [out] aParticles pointer to the list of data for all particles generated in current event
 * @param [in] aEventType enum that specifies which CALM option will be performed
 * 
 * 
 * 
 * @fn void CALM::SaveAllParticles_GLOBAL_REGGAE(int Nsum, vector4 *avec, double **XYZrand, ParticleDB *aPartDB, list<Particle> *aParticles)
 * @brief Pushes all particles into the array which will be saved to the file (for REGGAE GLOBAL)
 * @param [in] Nsum amount of all particles
 * @param [in] avec pointer to array of 4D vectors holding energy for each partible
 * @param [in] XYZrand array holding XYZ vertex coordinates got from CALM::GetXYZ
 * @param [in] aPartDB pointer to particle data base
 * @param [out] aParticles pointer to the list of data for all particles generated in current event
 * 
 * 
 * 
 * @fn void CALM::SaveAllParticles_MINIJETS_REGGAE(vector<double> *masses, vector<string> *names, vector4 *avec0, vector4 *avec1, double TotEnergy,
 *  double **XYZrand, ParticleDB *aPartDB, list<Particle> *aParticles)
 * @brief Pushes all particles into the array which will be saved to the file (for REGGAE MINIJETS)
 * @param [in] masses pointer to the array of vectors which will be filled by masses of particles for each jet
 * @param [in] names pointer to the array of vectors which will be filled by names of particles for each jet
 * @param [in] avec0 pointer to array of 4D vectors holding energy for each partible inside first jet
 * @param [in] avec1 pointer to array of 4D vectors holding energy for each partible inside second jet
 * @param [in] TotEnergy variable holding total energy
 * @param [in] XYZrand array holding XYZ vertex coordinates got from CALM::GetXYZ
 * @param [in] aPartDB pointer to particle data base
 * @param [out] aParticles pointer to the list of data for all particles generated in current event
 *  
 */
