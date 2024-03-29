/*This code is based on the code from THERMINATOR 2 generator, see details below*/

/********************************************************************************
 *                                                                              *
 *             THERMINATOR 2: THERMal heavy-IoN generATOR 2                     *
 *                                                                              *
 * Version:                                                                     *
 *      Release, 2.0.3, 1 February 2011                                         *
 *                                                                              *
 * Authors:                                                                     *
 *      Mikolaj Chojnacki   (Mikolaj.Chojnacki@ifj.edu.pl)                      *
 *      Adam Kisiel         (kisiel@if.pw.edu.pl)                               *
 *      Wojciech Broniowski (Wojciech.Broniowski@ifj.edu.pl)                    *
 *      Wojciech Florkowski (Wojciech.Florkowski@ifj.edu.pl)                    *
 *                                                                              *
 * Project homepage:                                                            *
 *      http://therminator2.ifj.edu.pl/                                         *
 *                                                                              *
 * For the detailed description of the program and further references           *
 * to the description of the model please refer to                              *
 * http://arxiv.org/abs/1102.0273                                               *
 *                                                                              *
 * This code can be freely used and redistributed. However if you decide to     *
 * make modifications to the code, please, inform the authors.                  *
 * Any publication of results obtained using this code must include the         *
 * reference to arXiv:1102.0273 and the published version of it, when           *
 * available.                                                                   *
 *                                                                              *
 ********************************************************************************/

#include <fstream>
#include <sstream>
#include <cstring>
#include "Configurator.h"
#include "THGlobal.h"
#include "Parser.h"

using namespace std;

Configurator::Configurator()
{
  mParameters.clear();
}

Configurator::Configurator(TString sConfigINI)
{
  mParameters.clear();

  _sConfigINI = sConfigINI;
}

Configurator::~Configurator()
{
  mParameters.clear();
}

void Configurator::ReadParameters()
{
  Parser* tParser;
  
  tParser = new Parser(_sConfigINI);
  tParser->ReadINI(this);
  delete tParser;
}

TString Configurator::GetParameter(const char* aKeyword)
{
  vector<Parameter>::iterator Iter;
  for (Iter = mParameters.begin(); Iter != mParameters.end(); Iter++)
    if (Iter->keyword == aKeyword) {
      PRINT_DEBUG_2("<Configurator::GetParameter>\tReturning value " << Iter->value << " for keyword " << Iter->keyword);
      return Iter->value;
    }
    return TString("");
}

void Configurator::AddParameter(Parameter* aPar)
{
  mParameters.push_back(*aPar);
}

int Configurator::PrintParameters()
{
  vector<Parameter>::iterator Iter;
  PRINT_MESSAGE("<Configurator::PrintParameters>");
  for (Iter = mParameters.begin(); Iter != mParameters.end(); Iter++)
    PRINT_MESSAGE("\tKeyword: " << Iter->keyword << "\tValue: " << Iter->value);
  return 0;
}
