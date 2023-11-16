#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TROOT.h"
#include "TRandom.h"
#include "particle.hpp"
#include "particle_type.hpp"
#include "resonance_type.hpp"

int Main() {
  TCanvas *canvas = new TCanvas("c", "canvas", 1280, 720);
  gRandom->SetSeed();

  Particle::AddParticleType("pione+", 0.13957, 1);     // pione+
  Particle::AddParticleType("pione-", 0.13957, -1);    // pione-
  Particle::AddParticleType("Kaone+", 0.49367, 1);     // Kaone+
  Particle::AddParticleType("Kaone-", 0.49367, -1);    // Kaone-
  Particle::AddParticleType("protone+", 0.93827, 1);   // protone+
  Particle::AddParticleType("protone-", 0.93827, -1);  // protone-
  Particle::AddParticleType("K*", 0.89166, 0, 0.050);  // risonanza K*

  std::vector<Particle> eventParticles{};
  for (int i = 0; i < 100; ++i) {
    Particle p;
    eventParticles.push_back(p);
  }

  // vettore che conterrà le "figlie" delle risonanze, poi lo appiccicheremo a
  // eventParticles
  std::vector<Particle> resonanceDaughters{};

  // istogrammi (numero bin e range da definire)
  TH1F *hParticleTypes = new TH1F("h1", "Particle Types", 7, 0, 7);
  TH1F *hPolar = new TH1F("h2", "Polar", 100, -100000, 100000);
  TH1F *hAzimuthal = new TH1F("h3", "Azimuthal", 100, -100000, 1000000);
  TH1F *hP = new TH1F("h4", "Impulse", 100, -100000, 10000000);
  TH1F *hTrsP = new TH1F("h5", "Trasverse impulse", 100, -100000, 100000);
  TH1F *hEnergy = new TH1F("h6", "Energy", 100, -1000000, 1000000);

  // istogrammi massa invariante
  // tutte le particelle
  TH1F *hInvariantMass0 = new TH1F("h7", "Invariant Mass", 100, -1000000,
                                   1000000);  // tutte le particelle
  hInvariantMass0->Sumw2();

  // particelle di segno discorde
  TH1F *hInvariantMass1 =
      new TH1F("h8", "Invariant Mass Opposite Sign", 100, -1000000, 1000000);
  hInvariantMass1->Sumw2();

  // particelle di segno concorde
  TH1F *hInvariantMass2 =
      new TH1F("h9", "Invariant Mass Same Sign", 100, -1000000, 1000000);
  hInvariantMass2->Sumw2();

  // pione+/Kaone- or pione-/Kaone+
  TH1F *hInvariantMass3 =
      new TH1F("h10", "Invariant Mass pione+/Kaone- or pione-/Kaone+", 100,
               -100000, 1000000);
  hInvariantMass3->Sumw2();

  // pione+/Kaone+ or pione-/Kaone-
  TH1F *hInvariantMass4 =
      new TH1F("h11", "Invariant Mass pione+/Kaone+ or pione-/Kaone-", 100,
               -1000000, 1000000);
  hInvariantMass4->Sumw2();

  // figlie delle risonanze
  TH1F *hInvariantMass5 =
      new TH1F("h12", "Invariant Mass Decay Particles", 100, -1000000, 1000000);
  hInvariantMass5->Sumw2();

  // inizio dei 10e5 eventi
  for (int i = 0; i < 10e5; ++i) {
    // generiamo le particelle e nient'altro, degli istogrammi ci preoccupiamo
    // in un secondo momento
    for (int j = 0; j < 100; ++j) {
      double phi = gRandom->Uniform(0, 2 * M_PI);
      hAzimuthal->Fill(phi);

      double theta = gRandom->Uniform(0, M_PI);
      hPolar->Fill(theta);

      double p = gRandom->Exp(1.);
      eventParticles[j].SetP(p * std::sin(theta) * std::cos(phi),
                             p * std::sin(theta) * std::sin(phi),
                             p * std::cos(theta));

      int rand = (gRandom->Uniform(1., 101.));
      if (rand <= 80) {
        if (rand <= 40) {
          eventParticles[j].SetIndex("pione+");
        } else {
          eventParticles[j].SetIndex("pione-");
        }
      } else if (rand <= 90) {
        if (rand <= 85) {
          eventParticles[j].SetIndex("Kaone+");
        } else {
          eventParticles[j].SetIndex("Kaone-");
        }
      } else if (rand <= 99) {
        if (std::rand() % 2) {
          eventParticles[j].SetIndex("protone+");
        } else {
          eventParticles[j].SetIndex("protone-");
        }
      } else if (rand == 100) {
        // non avevamo capito la logica di Decay2body()
        eventParticles[j].SetIndex("K*");

        Particle p1{};
        Particle p2{};

        if (std::rand() % 2) {
          p1.SetIndex("pione+");
          p2.SetIndex("Kaone-");
        } else {
          p1.SetIndex("pione-");
          p2.SetIndex("Kaone+");
        }

        eventParticles[j].Decay2body(p1, p2);

        resonanceDaughters.push_back(p1);
        resonanceDaughters.push_back(p2);
      }
    }  // fine della generazione delle ~100 particelle

    // inserisco le "figlie" delle risonanze in coda a eventParticles
    eventParticles.insert(eventParticles.end(), resonanceDaughters.begin(),
                          resonanceDaughters.end());

    // riempimento istogrammi
    for (Particle particle : eventParticles) {
      hParticleTypes->Fill(particle.GetIndex());
      hP->Fill(particle.GetP());
      hTrsP->Fill(particle.GetTrsP());
      hEnergy->Fill(particle.TotalEnergy());
    }

    // riempimento istogrammi massa invariante
    for (int i = 0; i < static_cast<int>(eventParticles.size()); ++i) {
      for (int j = 0; static_cast<int>(eventParticles.size()); ++j) {
        // tutte le particelle
        hInvariantMass0->Fill(
            eventParticles[i].InvariantMass(eventParticles[j]));

        int indexI = eventParticles[i].GetIndex();
        int indexJ = eventParticles[j].GetIndex();

        // escludiamo le risonanze
        if (indexI != 6 && indexJ != 6) {
          // segno discorde
          if (indexI % 2 != indexJ % 2) {
            hInvariantMass1->Fill(
                eventParticles[i].InvariantMass(eventParticles[j]));
          }
          // segno concorde
          else {
            hInvariantMass2->Fill(
                eventParticles[i].InvariantMass(eventParticles[j]));
          }

          // pione+/Kaone- e pione-/Kaone+ (i.e. pioni e Kaoni con segni
          // discordi)
          if (((indexI == 0 || indexJ == 0) && (indexI + indexJ == 3)) ||
              ((indexI == 1 || indexJ == 1) && (indexI + indexJ == 2))) {
            hInvariantMass4->Fill(
                eventParticles[i].InvariantMass(eventParticles[j]));
          }

          // pione+/Kaone+ e pione-/Kaone- (i.e. pioni e Kaoni con segni
          // concordi)
          if (((indexI == 0 || indexJ == 0) && (indexI + indexJ == 2)) ||
              ((indexI == 1 || indexJ == 1) && (indexI + indexJ == 4))) {
            hInvariantMass5->Fill(
                eventParticles[i].InvariantMass(eventParticles[j]));
          }
        }
      }
    }

    // riempimento istogramma  massa invariante delle figlie della risonanza
    // (solo figlie provenienti dalla stessa madre, per questo i += 2)
    for (int i = 0; i < static_cast<int>(resonanceDaughters.size()); i += 2) {
      hInvariantMass5->Fill(
          resonanceDaughters[i].InvariantMass(resonanceDaughters[i + 1]));
    }
  }  // fine dei 10e5 eventi

  // inserimento degli istogrammi in un file ROOT - per ora lasciamo stare
  /* TFile *file = new TFile("histograms.root");

  hParticleTypes->Write();
  hPolar->Write();
  hAzimuthal->Write();
  hP->Write();
  hTrsP->Write();
  hEnergy->Write();
  hInvariantMass0->Write();
  hInvariantMass1->Write();
  hInvariantMass2->Write();
  hInvariantMass3->Write();
  hInvariantMass4->Write();
  hInvariantMass5->Write();

  file->Close(); */

  // per testare il funzionamento del codice
  hParticleTypes->Draw();
  hInvariantMass5->Draw();

  return 0;
}
