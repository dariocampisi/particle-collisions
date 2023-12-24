#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TRandom.h"
#include "particle.hpp"
#include "particle_type.hpp"
#include "resonance_type.hpp"

void HistoAnalysis() {
  TFile *file = new TFile("histograms.root");

  TH1F *hParticleTypes = (TH1F *)file->Get("A");
  TH1F *hPolar = (TH1F *)file->Get("B");
  TH1F *hAzimuthal = (TH1F *)file->Get("C");
  TH1F *hImpulse = (TH1F *)file->Get("D");
  TH1F *hTrsImpulse = (TH1F *)file->Get("E");
  TH1F *hEnergy = (TH1F *)file->Get("F");
  TH1F *hIMAll = (TH1F *)file->Get("G");
  TH1F *hIMOppositeCharge = (TH1F *)file->Get("H");
  TH1F *hIMSameCharge = (TH1F *)file->Get("I");
  TH1F *hIMOppositePK = (TH1F *)file->Get("J");
  TH1F *hIMSamePK = (TH1F *)file->Get("K");
  TH1F *hIMDecayParticles = (TH1F *)file->Get("L");

  std::vector<TH1F *> histoVector{
      hParticleTypes, hPolar,        hAzimuthal, hImpulse,
      hTrsImpulse,    hEnergy,       hIMAll,     hIMOppositeCharge,
      hIMSameCharge,  hIMOppositePK, hIMSamePK,  hIMDecayParticles};

  // NUMERO DI INGRESSI PER OGNI ISTOGRAMMA
  std::cout << "\n\n\n********** ENTRIES **********\n";
  for (auto histo : histoVector) {
    std::cout << "\n " << histo->GetTitle() << ": " << histo->GetEntries();
  }

  // DISTRIBUZIONE DEI TIPI DI PARTICELLE
  std::cout << "\n\n********** PARTICLE TYPES **********\n\n";

  std::cout << " Pions (+): " << hParticleTypes->GetBinContent(1) << " i.e. "
            << 100 * hParticleTypes->GetBinContent(1) /
                   hParticleTypes->GetEntries()
            << "%\n";
  std::cout << " Pions (-): " << hParticleTypes->GetBinContent(2) << " i.e. "
            << 100 * hParticleTypes->GetBinContent(2) /
                   hParticleTypes->GetEntries()
            << "%\n";
  std::cout << " Kaons (+): " << hParticleTypes->GetBinContent(3) << " i.e. "
            << 100 * hParticleTypes->GetBinContent(3) /
                   hParticleTypes->GetEntries()
            << "%\n";
  std::cout << " Kaons (-): " << hParticleTypes->GetBinContent(4) << " i.e. "
            << 100 * hParticleTypes->GetBinContent(4) /
                   hParticleTypes->GetEntries()
            << "%\n";
  std::cout << " Protons (+): " << hParticleTypes->GetBinContent(5) << " i.e. "
            << 100 * hParticleTypes->GetBinContent(5) /
                   hParticleTypes->GetEntries()
            << "%\n";
  std::cout << " Protons (-): " << hParticleTypes->GetBinContent(6) << " i.e. "
            << 100 * hParticleTypes->GetBinContent(6) /
                   hParticleTypes->GetEntries()
            << "%\n";
  std::cout << " Resonances: " << hParticleTypes->GetBinContent(7) << " i.e. "
            << 100 * hParticleTypes->GetBinContent(7) /
                   hParticleTypes->GetEntries()
            << "%\n";

  // FITTING
  // angolo polare
  TF1 *polarFit = new TF1("polar fit function", "[0]", 0., M_PI);
  hPolar->Fit(polarFit, "Q0");

  // angolo azimutale
  TF1 *azimuthalFit = new TF1("azimuthal fit function", "[0]", 0., 2 * M_PI);
  hAzimuthal->Fit(azimuthalFit, "Q0");

  // impulso
  TF1 *impulseFit = new TF1("impulse fit function", "expo", 0., 7.);
  hImpulse->Fit(impulseFit, "Q0");

  std::cout << "\n********** FITTING **********\n\n";

  std::cout << " Polar Angle\n\n"
            << "  Parameter: " << polarFit->GetParameter(0) << " ± "
            << polarFit->GetParError(0) << "\n  Reduced Chi Square: "
            << polarFit->GetChisquare() / polarFit->GetNDF()
            << "\n  Probability: " << polarFit->GetProb();

  std::cout << "\n\n Azimuthal Angle\n\n"
            << "  Parameter: " << azimuthalFit->GetParameter(0) << " ± "
            << azimuthalFit->GetParError(0) << "\n  Reduced Chi Square: "
            << azimuthalFit->GetChisquare() / azimuthalFit->GetNDF()
            << "\n  Probability: " << azimuthalFit->GetProb();

  std::cout << "\n\n Impulse\n\n"
            << "  Width: " << impulseFit->GetParameter(0) << " ± "
            << impulseFit->GetParError(0)
            << "\n  Mean: " << impulseFit->GetParameter(1) << " ± "
            << impulseFit->GetParError(1) << "\n  Reduced Chi Square: "
            << impulseFit->GetChisquare() / impulseFit->GetNDF()
            << "\n  Probability: " << impulseFit->GetProb();

  // SOTTRAZIONE ISTOGRAMMI
  TH1F *hSubtraction1 = new TH1F(*hIMOppositeCharge);
  hSubtraction1->Add(hIMOppositeCharge, hIMSameCharge, 1., -1.);
  TF1 *sub1Fit = new TF1("first subtraction fit", "gaus", 0., 8.);
  hSubtraction1->Fit(sub1Fit, "Q0");

  TH1F *hSubtraction2 = new TH1F(*hIMOppositePK);
  hSubtraction2->Add(hIMOppositePK, hIMSamePK, 1., -1.);
  TF1 *sub2Fit = new TF1("second subtraction fit", "gaus", 0., 8.);
  hSubtraction2->Fit(sub2Fit, "Q0");

  std::cout << "\n\n********** OPERATIONS ON HISTOGRAMS **********\n\n";

  std::cout << " Inv. Mass - opposite charge minus Inv. Mass - same charge\n\n"
            << "  Width: " << sub1Fit->GetParameter(0) << " ± "
            << sub1Fit->GetParError(0)
            << "\n  Mean (i.e. resonance mass): " << sub1Fit->GetParameter(1)
            << " ± " << sub1Fit->GetParError(1)
            << "\n  Std Dev (i.e. resonance width): "
            << sub1Fit->GetParameter(2) << " ± " << sub1Fit->GetParError(2)
            << "\n  Reduced Chi Square: "
            << sub1Fit->GetChisquare() / sub1Fit->GetNDF()
            << "\n  Fit Probability: " << sub1Fit->GetProb();

  std::cout << "\n\n Inv. Mass - opposite charge pions and kaons minus Inv. "
               "Mass - same charge pions and kaons\n\n"
            << "  Width: " << sub2Fit->GetParameter(0) << " ± "
            << sub2Fit->GetParError(0)
            << "\n  Mean (i.e. resonance mass): " << sub2Fit->GetParameter(1)
            << " ± " << sub2Fit->GetParError(1)
            << "\n  Std Dev (i.e. resonance width): "
            << sub2Fit->GetParameter(2) << " ± " << sub2Fit->GetParError(2)
            << "\n  Reduced Chi Square: "
            << sub2Fit->GetChisquare() / sub2Fit->GetNDF()
            << "\n  Fit Probability: " << sub2Fit->GetProb() << "\n\n\n";
}



// nel seguito solo porcherie scritte da giovanni brandi da cui mi dissocio completamente e senza riserve

/*TH1F *h1 = (TH1F *)file->Get("h1");
int expectedEntries = 10E7;
int entries = histo1->GetEntries();
double ratio = entries / expectedEntries;
std::cout << "\n\n******** Analisi numero di ingressi ********";
std::cout << "\nExpected entries: " << expectedEntries;
std::cout << "\nEntries: " << entries;
std::cout << "\nRatio:" << histo1->GetNbinsX();

if (ratio > 1.1) {
  std::cout << "\nNumero di ingressi troppo alto";
} else if (ratio < 0.9) {
  std::cout << "\nNumero di ingressi troppo basso";
} else {
  std::cout << "\nNumero di ingresi atteso";
}

std::cout << "\n******** Analisi percentuali di particelle ********";

TH1F *hParticleTypes = (TH1F *)file->Get("hParticleTypes");
int expectedPercentages[7] = {
    0.4, 0.4, 0.05, 0.05, 0, 0, 0};  // allora levo anche getnbins sotto
for (int i = 0; i < hParticleTypes->GetNbinsX(); i++) {
  double upperBound =
      expectedPercentages[i] + hParticleTypes->GetBinError(i) / entries;
  double lowerBound =
      expectedPercentages[i] - hParticleTypes->GetBinError(i) / entries;
  double percentage = hParticleTypes->GetBinContent(i) / entries;
  std::cout << "\nTipo di particella: " << i;
  std::cout << "\nPercentuale sul totale: " << percentage << " +- "
            << hParticleTypes->GetBinError(i) / entries;
  if (percentage > upperBound) {
    std::cout << "\nLa percentuale è maggiore di quella attesa.";
  } else if (percentage < lowerBound) {
    std::cout << "\nLa percentuale è minore di quella attesa.";
  } else {
    std::cout << "\nLa percentuale è coerente con quella attesa.";
  }
}*/