#include "Pythia8/Pythia.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"

using namespace Pythia8;

void saveHistogramToFile(TH1F* hist, const char* filename) {
    TFile outputFile(filename, "RECREATE");
    hist->Write();
    outputFile.Close();
}

int main() {
    Pythia pythia;

    // Ustawienie parametrów symulacji.
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Main:numberOfEvents = 100000");
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 7000");
    pythia.readString("SoftQCD:nonDiffractive = on");
    pythia.readString("SoftQCD:singleDiffractive = on");
    pythia.readString("SoftQCD:doubleDiffractive = on");
    pythia.readString("PhaseSpace:pTHatMin = 0");
    pythia.readString("PhaseSpace:pTHatMax = 7000");
    pythia.readString("PhaseSpace:mHatMin = 0");
    pythia.readString("PhaseSpace:mHatMax = 7000");
    pythia.readString("ParticleDecays:limitTau0 = On");
    pythia.readString("ParticleDecays:tau0Max = 10.0");
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 0");

    // Inicjalizacja.
    pythia.init();

    // Histogram do przechowywania rozkładów. 
    TH1F *energyDistribution = new TH1F("Energy distribution", "Energy distribution", 9000, 0, 1100);
    TH1F *lambdasDistribution= new TH1F("Lambdas distribution", "Lambdas distribution", 10, 0, 9);
    TH1F *nucleonsDistribution = new TH1F("Nucleons distribution", "Nucleons distribution", 28, 0, 27);
    TH1F *kaonsDistribution = new TH1F("Kaons distribution", "Kaons distribution", 38, 0, 37);
    TH1F *pionsDistribution = new TH1F("Pions distribution", "Pions distribution", 6000, 0, 120);

    // Generowanie zderzeń.
    for (int iEvent = 0; iEvent < 100000; ++iEvent) {
        if (!pythia.next()) continue;
        
        int lambdaCount = 0;
        int nucleonCount = 0;
        int kaonCount = 0;
        int pionCount = 0;

        // Przeszukiwanie wszystkich cząstek w zderzeniu.
        for (int i = 0; i < pythia.event.size(); ++i) {
            int pdgId = pythia.event[i].id();
            double rapidity = pythia.event[i].eta();
            double energy = pythia.event[i].e();
            int status = pythia.event[i].status();
                
            if (status >= 81 && status <=89) {
                if (pdgId == 3122 || pdgId == -3122)
                {
                	lambdaCount++;
                	energyDistribution->Fill(energy);
                }
            	else if (pdgId == 311 || pdgId == -311 || pdgId == 321 || pdgId == -321)
            	{
            		kaonCount++;
            		energyDistribution->Fill(energy);
            	} 
            	else if (pdgId == 211 || pdgId == -211)
            	{
            		pionCount++;
            		energyDistribution->Fill(energy);
            	}  
            	else if (pdgId == 2212 || pdgId == -2212 || pdgId == 2112 || pdgId == -2112)
            	{
            		nucleonCount++;
            		energyDistribution->Fill(energy);
            	}
            }
        }
        
        lambdasDistribution->Fill(lambdaCount);
        nucleonsDistribution->Fill(nucleonCount);
        kaonsDistribution->Fill(kaonCount);
        pionsDistribution->Fill(pionCount);
    }
    
    pionsDistribution->Scale(1.0 / pionsDistribution->Integral(), "width");
    kaonsDistribution->Scale(1.0 / kaonsDistribution->Integral(), "width");
    nucleonsDistribution->Scale(1.0 / nucleonsDistribution->Integral(), "width");  
    lambdasDistribution->Scale(1.0 / lambdasDistribution->Integral(), "width");
    energyDistribution->Scale(1.0 / energyDistribution->Integral(), "width");

    saveHistogramToFile(pionsDistribution, "pionsDistribution.root");
    saveHistogramToFile(kaonsDistribution, "kaonsDistribution.root");
    saveHistogramToFile(nucleonsDistribution, "nucleonsDistribution.root");
    saveHistogramToFile(lambdasDistribution, "lambdasDistribution.root");
    saveHistogramToFile(energyDistribution, "energyDistribution.root");
    
    TFile outputFile("multiplicity_histograms.root", "RECREATE");
        
    pionsDistribution->Write();
    kaonsDistribution->Write();
    nucleonsDistribution->Write();
    lambdasDistribution->Write();
    energyDistribution->Write();

    outputFile.Close();

    TCanvas c("c", "Particle Multiplicity Distributions", 800, 600);
    gStyle->SetOptStat(000);
    c.SetLogy(1);
    
    pionsDistribution->GetXaxis()->SetTitle("N_{ch}");
    pionsDistribution->GetYaxis()->SetTitle("dN/dN_{ch}");
    pionsDistribution->Draw("P");
    c.SaveAs("pionsDistribution.png");

    kaonsDistribution->GetXaxis()->SetTitle("N_{ch}");
    kaonsDistribution->GetYaxis()->SetTitle("dN/dN_{ch}");
    kaonsDistribution->Draw("P");
    c.SaveAs("kaonsDistribution.png");

    nucleonsDistribution->GetXaxis()->SetTitle("N_{ch}");
    nucleonsDistribution->GetYaxis()->SetTitle("dN/dN_{ch}");
    nucleonsDistribution->Draw("P");
    c.SaveAs("nucleonsDistribution.png");

    lambdasDistribution->GetXaxis()->SetTitle("N_{ch}");
    lambdasDistribution->GetYaxis()->SetTitle("dN/dN_{ch}");
    lambdasDistribution->Draw("P");
    c.SaveAs("lambdasDistribution.png");
    
    energyDistribution->GetXaxis()->SetTitle("E (GeV)");
    energyDistribution->GetYaxis()->SetTitle("dN/dN_{ch}");
    energyDistribution->Draw("P");
    c.SaveAs("energyDistribution.png");
      
    pythia.stat();

    return 0;
}
