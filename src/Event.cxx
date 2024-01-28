/*This code is based on the code from THERMINATOR 2 generator, see details below*/
/********************************************************************************
 *																				*
 *						 THERMINATOR 2: THERMal heavy-IoN generATOR 2			*
 *																				*
 * Version:																		*
 *			Release, 2.0.3, 1 February 2011										*
 *																				*
 * Authors:																		*
 *			Mikolaj Chojnacki	 (Mikolaj.Chojnacki@ifj.edu.pl)					*
 *			Adam Kisiel				 (kisiel@if.pw.edu.pl)						*
 *			Wojciech Broniowski (Wojciech.Broniowski@ifj.edu.pl)				*
 *			Wojciech Florkowski (Wojciech.Florkowski@ifj.edu.pl)				*
 *																				*
 * Project homepage:															*
 *			http://therminator2.ifj.edu.pl/										*
 *																				*
 * For the detailed description of the program and further references			*
 * to the description of the model please refer to								*
 * http://arxiv.org/abs/1102.0273												*
 *																				*
 * This code can be freely used and redistributed. However if you decide to		*
 * make modifications to the code, please, inform the authors.					*
 * Any publication of results obtained using this code must include the			*
 * reference to arXiv:1102.0273 and the published version of it, when			*
 * available.																	*
 *																				*
 ********************************************************************************/

#include <sstream>
#include <TDatime.h>
#include "Crc32.h"
#include "Configurator.h"
#include "Event.h"
#include "THGlobal.h"

extern Configurator* sMainConfig;
extern TString	sTimeStamp;
extern int	sRandomize;

using namespace std;

Event::Event()
   : mPartDB(0), mCALM(0), mRandom(0), mDistribution(0), mMultMin(10), mMultMax(20), mEnergy(30), mEventType(GLOBAL), mJetNumber(2)
{
	mMultiplicities.clear();
	Reset();
}

Event::Event(ParticleDB* aDB, CALM* aCALM)
: mPartDB(aDB), 
 mCALM(aCALM), mDistribution(0), mMultMin(10), mMultMax(20), mEnergy(30)
{ 
	mRandom = new TRandom2();
#ifdef _ROOT_4_
	mRandom->SetSeed2(31851, 14327);
#else
	mRandom->SetSeed(31851);
#endif
	mMultiplicities.clear();
	mMultiplicities.resize(mPartDB->GetParticleTypeCount());
	Reset();
	ReadParameters();
}

Event::~Event()
{
	mParticles.clear();
	mMultiplicities.clear();
	delete mRandom;
}

void Event::Reset(int aEventIter)
{
	ostringstream oss;
	Crc32 tEventID;
	
	mParticles.clear();
	Particle::ZeroEID();
	
	oss << sTimeStamp.Data() << "Event: " << aEventIter;
	tEventID.Update(oss.str().data(), oss.str().length());
	tEventID.Finish(); 
	mEventID = tEventID.GetValue();
}

list<Particle>* Event::GetParticleList()
{
	return &mParticles;
}

ParticleDB* Event::GetParticleDB() const
{
	return mPartDB;
}

unsigned int Event::GetEventID() const
{
	return mEventID;
}

void Event::GeneratePrimordials(int aSeed)
{ 
#ifdef _ROOT_4_
	if (aSeed) mRandom->SetSeed2(aSeed, (aSeed*2) % (7*11*23*31));
#else
	if (aSeed) mRandom->SetSeed(aSeed);
#endif
	int control=99;
	//GenerateMultiplicities();
	do
	  {
	    control = mCALM->GenerateParticles(mPartDB, mMultMin, mMultMax, mEnergy, &mParticles, mEventType, mJetNumber);
	  }
	while( control == 99 );
		
}

void Event::GenerateMultiplicities()
{
	if(mDistribution == 0) { // Poisson
		for (int tIter=0; tIter<mPartDB->GetParticleTypeCount(); tIter++)
			mMultiplicities[tIter] = mRandom->Poisson(mPartDB->GetParticleType(tIter)->GetMultiplicity());
	} else if(mDistribution == 1) { // Negative Binomial
		for (int tIter=0; tIter<mPartDB->GetParticleTypeCount(); tIter++)
			mMultiplicities[tIter] = 0; // HOW?
	}
}

void Event::Randomize()
{
	TDatime tDate;

#ifdef _ROOT_4_
	mRandom->SetSeed2(tDate.Get() / 2 * 3, tDate.Get() / 11 * 9);
#else
	mRandom->SetSeed(tDate.Get() / 2 * 3);
#endif
}

void Event::ReadParameters()
{
	TString tDistribution; 
	std::string tMultMin, tMultMax, tEnergy, tEventType, tJetNumber;
	int tMultMinInt, tMultMaxInt, tJetNumberInt;
   int tEventTypeEnum;
	double tEnergyDouble;
	std::stringstream tConvert;
	try {
		tDistribution	= sMainConfig->GetParameter("MultiplicityDistribution");
		if (tDistribution.Contains("NegativeBinomial"))
			mDistribution = 1;
		tMultMin	= sMainConfig->GetParameter("MultiplicityMin");
		tMultMax	= sMainConfig->GetParameter("MultiplicityMax");
		tEnergy	= sMainConfig->GetParameter("GenbodEnergy");
		tEventType	= sMainConfig->GetParameter("EventType");
		tJetNumber = sMainConfig->GetParameter("NumberOfJets");
		tConvert<<tMultMin<<' '<<tMultMax<<' '<<tEnergy<<' '<<tEventType<<' '<<tJetNumber;
		tConvert>>tMultMinInt>>tMultMaxInt>>tEnergyDouble>>tEventTypeEnum>>tJetNumberInt;
		if(tMultMinInt<tMultMaxInt)
		{
			mMultMin = tMultMinInt;
			mMultMax = tMultMaxInt;
		}
      mEventType = (eEventType) tEventTypeEnum;
		if(tEnergyDouble>0)	mEnergy = tEnergyDouble;
		if(tJetNumberInt>1 && tJetNumberInt<5) mJetNumber = tJetNumberInt;
		PRINT_DEBUG_1("<Event::ReadParameters>\tSetting multiplicity range ("<<mMultMin<<","<<mMultMax<<") and energy "<<mEnergy);
	}
	catch (TString tError) {
		PRINT_DEBUG_1("<Event::ReadParameters>\tUsing default multiplicity distribution: Poissonian");
	}
}
