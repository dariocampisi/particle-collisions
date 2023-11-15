#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TRandom.h"
#include "particle.hpp"
#include "particle_type.hpp"
#include "resonance_type.hpp"

void ReadHistograms() {
  TFile *file = new TFile("histograms.root");

  for (int i = 1; i <= 12; ++i) {
    char c = static_cast<char>(i + '0');
    char *pt = &c;
    TH1F *h = (TH1F *)file->Get('h' + pt);
  }

  TH1F *h1 = (TH1F *)file->Get("h1");
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
  }

  std::cout << "\n**** Verifica consistenza distribuzione direzioni angolari "
               "con dist. uniforme*****";
  TF1 *polarFit = new TF1("polarFit", "[0]", 0, M_PI);
  TH1F *hPolar = (TH1F *)file->Get("hPolar");
  hPolar->Fit(polarFit);
  TF1 *polarFitFunc = hPolar->GetFunction("polarFit");
  std::cout << "\nComponente polare:";
  std::cout << "\n    Parametro del fit: " << polarFitFunc->GetParameter(0);
  std::cout << "\n    Chi quadro: " << polarFitFunc->GetChisquare();
  std::cout << "\n    Gradi di libertà: " << polarFitFunc->GetNDF();
  std::cout << "\n    Probabilità di ottenere il fit: ";

  TF1 *azimuthalFit = new TF1("f1", "[0]", 0, 2 * M_PI);
  TH1F *hAzimuthal = (TH1F *)file->Get("hAzimuthal");
  hAzimuthal->Fit(azimuthalFit);
  TF1 *azimuthalFitFunc = hAzimuthal->GetFunction("azimuthalFit");
  std::cout << "\n\nComponente azimuthale:";
  std::cout << "\n    Parametro del fit: " << azimuthalFitFunc->GetParameter(0);
  std::cout << "\n    Chi quadro: " << azimuthalFitFunc->GetChisquare();
  std::cout << "\n    Gradi di libertà: " << azimuthalFitFunc->GetNDF();
  std::cout << "\n    Probabilità di ottenere il fit: ";
}
