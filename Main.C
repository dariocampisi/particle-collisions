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
  gRandom->SetSeed(23423);

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

  std::vector<Particle> resonanceDaughters{};

  TH1F *hParticleTypes = new TH1F("h1", "Particle Types", 7, 0, 7);
  TH1F *hPolar = new TH1F("h2", "Polar", 100, -100000, 100000);
  TH1F *hAzimuthal = new TH1F("h3", "Azimuthal", 100, -100000, 1000000);
  TH1F *hP = new TH1F("h4", "Impulse", 100, -100000, 10000000);
  TH1F *hTrsP = new TH1F("h5", "Trasverse impulse", 100, -100000, 100000);
  TH1F *hEnergy = new TH1F("h6", "Energy", 100, -1000000, 1000000);

  TH1F *hInvariantMass0 = new TH1F("h7", "Invariant Mass", 100, -1000000,
                                   1000000);  // tutte le particelle
  TH1F *hInvariantMass1 =
      new TH1F("h8", "Invariant Mass Opposite Sign", 100, -1000000,
               1000000);  // segno discorde
  TH1F *hInvariantMass2 =
      new TH1F("h9", "Invariant Mass Same Sign", 100, -1000000,
               1000000);  // segno concorde
  TH1F *hInvariantMass3 = new TH1F(
      "h10", "Invariant Mass pione+/Kaone- or pione-/Kaone+", 100, -100000,
      1000000);  // pione+/Kaone- or pione-/Kaone+
  TH1F *hInvariantMass4 = new TH1F(
      "h11", "Invariant Mass pione+/Kaone+ or pione-/Kaone-", 100, -1000000,
      1000000);  // pione+/Kaone+ or pione-/Kaone-
  TH1F *hInvariantMass5 =
      new TH1F("h12", "Invariant Mass Decay Particles", 100, -1000000,
               1000000);  // figlie delle risonanze

  hInvariantMass0->Sumw2();
  hInvariantMass1->Sumw2();
  hInvariantMass2->Sumw2();
  hInvariantMass3->Sumw2();
  hInvariantMass4->Sumw2();
  hInvariantMass5->Sumw2();

  for (int i = 0; i < 10e5; ++i) {
    for (int j = 0; j < 100; ++j) {
      double phi = gRandom->Uniform(0, 2 * M_PI);
      hAzimuthal->Fill(phi);
      double theta = gRandom->Uniform(0, M_PI);
      hPolar->Fill(theta);

      // double p = gRandom->Exp(1.);
      // hP->Fill(p);
      // double px = p * sin(theta) * cos(phi);
      // double py = p * sin(theta) * sin(phi);
      // double pz = p * cos(theta);
      // eventParticles[j].SetP(px, py, pz);

      // hTrsP->Fill(std::sqrt(px * px + py * py));

      // hEnergy->Fill(eventParticles[j].TotalEnergy());

      int rand = static_cast<int>(gRandom->Uniform(1., 101.));
      if (rand <= 80) {
        if (rand <= 40) {
          eventParticles[j].SetIndex("pione+");
          hParticleTypes->Fill(eventParticles[j].GetIndex());
        } else {
          eventParticles[j].SetIndex("pione-");
          hParticleTypes->Fill(eventParticles[j].GetIndex());
        }
      } else if (rand <= 90) {
        if (rand <= 85) {
          eventParticles[j].SetIndex("Kaone+");
          hParticleTypes->Fill(eventParticles[j].GetIndex());
        } else {
          eventParticles[j].SetIndex("Kaone-");
          hParticleTypes->Fill(eventParticles[j].GetIndex());
        }
      } else if (rand <= 99) {
        if (std::rand() % 2) {
          eventParticles[j].SetIndex("protone+");
          hParticleTypes->Fill(eventParticles[j].GetIndex());
        } else {
          eventParticles[j].SetIndex("protone-");
          hParticleTypes->Fill(eventParticles[j].GetIndex());
        }
      } else if (rand == 100) {
        Printf("ECCOLO");

        eventParticles[j].SetIndex("K*");
        hParticleTypes->Fill(eventParticles[j].GetIndex());

        Particle p1{};
        Particle p2{};
        eventParticles[j].Decay2body(p1, p2);

        if (std::rand() % 2) {
          p1.SetIndex("pione+");
          p2.SetIndex("Kaone-");
          resonanceDaughters.push_back(p1);
          resonanceDaughters.push_back(p2);
        } else {
          p1.SetIndex("pione-");
          p2.SetIndex("Kaone+");
          resonanceDaughters.push_back(p1);
          resonanceDaughters.push_back(p2);
        }
      }
    }  // fine della generazione delle ~100 particelle

    // inserisco le "figlie" della risonanza in coda al vettore delle particelle
    eventParticles.insert(eventParticles.end(), resonanceDaughters.begin(),
                          resonanceDaughters.end());

    // riempimento istogrammi massa invariante
    for (int i = 0; i < static_cast<int>(eventParticles.size()); ++i) {
      for (int j = i + 1; j < static_cast<int>(eventParticles.size()); ++j) {
        hInvariantMass0->Fill(
            eventParticles[i].InvariantMass(eventParticles[j]));

        int index1 = eventParticles[i].GetIndex();
        int index2 = eventParticles[j].GetIndex();

        // escludiamo le risonanze
        if (index1 != 6 && index2 != 6) {
          if (index1 % 2 != index2 % 2) {
            hInvariantMass1->Fill(eventParticles[i].InvariantMass(
                eventParticles[j]));  // segni concordi
          } else {
            hInvariantMass2->Fill(eventParticles[i].InvariantMass(
                eventParticles[j]));  // segni discordi
          }

          if (index1 + index2 == 3) {
            if (index1 % 2 != index2 % 2) {
              hInvariantMass4->Fill(eventParticles[i].InvariantMass(
                  eventParticles[j]));  // pione+/Kaone- e pione-/Kaone+
            } else {
              hInvariantMass5->Fill(eventParticles[i].InvariantMass(
                  eventParticles[j]));  // pione+/Kaone+ e pione-/Kaone-
            }
          }
        }
      }
    }
    // massa invariante delle figlie delle risonanze (solo figlie provenienti
    // dalla stessa madre)
    for (int i = 0; i < static_cast<int>(resonanceDaughters.size()); i += 2) {
      hInvariantMass5->Fill(
          resonanceDaughters[i].InvariantMass(resonanceDaughters[i + 1]));
    }
  }  // fine dei 10e5 eventi

  // inserimento degli istogrammi in un file ROOT
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

  return 0;
}
