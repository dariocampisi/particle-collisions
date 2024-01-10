#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TRandom.h"
#include "particle.hpp"
#include "particle_type.hpp"
#include "resonance_type.hpp"

void HistoAnalysis() {
  TFile *file = new TFile("histograms.root");

  TH1F *hParticleTypes = (TH1F *)file->Get("Particle Types");
  TH1F *hPolar = (TH1F *)file->Get("Polar Angle");
  TH1F *hAzimuthal = (TH1F *)file->Get("Azimuthal Angle");
  TH1F *hImpulse = (TH1F *)file->Get("Impulse");
  TH1F *hTrsImpulse = (TH1F *)file->Get("Trasverse Impulse");
  TH1F *hEnergy = (TH1F *)file->Get("Energy");
  TH1F *hIMAll = (TH1F *)file->Get("Invariant Mass");
  TH1F *hIMOppositeCharge = (TH1F *)file->Get("Inv. Mass - opposite charge");
  TH1F *hIMSameCharge = (TH1F *)file->Get("Inv. Mass - same charge");
  TH1F *hIMOppositePK =
      (TH1F *)file->Get("Inv. Mass - opposite charge pions and kaons");
  TH1F *hIMSamePK =
      (TH1F *)file->Get("Inv. Mass - same charge pions and kaons");
  TH1F *hIMDecayParticles = (TH1F *)file->Get("Inv. Mass - decay particles");

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
