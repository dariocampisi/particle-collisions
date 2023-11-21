#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TROOT.h"
#include "TRandom.h"
#include "particle.hpp"
#include "particle_type.hpp"
#include "resonance_type.hpp"

int Main() {
  gRandom->SetSeed(314234);

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

  std::vector<Particle> decayParticles{};

  // istogrammi (numero bin e range da definire)
  TH1F *hParticleTypes = new TH1F("A", "Particle Types", 7, 0., 7.);
  TH1F *hPolar = new TH1F("B", "Polar Angle", 180, 0., M_PI);
  TH1F *hAzimuthal = new TH1F("C", "Azimuthal Angle", 360, 0., 2 * M_PI);
  TH1F *hImpulse = new TH1F("D", "Impulse", 500, 0., 7.);
  TH1F *hTrsImpulse = new TH1F("E", "Trasverse Impulse", 500, 0., 7.);
  TH1F *hEnergy = new TH1F("F", "Energy", 100, -10, 10);

  // istogrammi massa invariante
  TH1F *hIMAll = new TH1F("G", "Invariant Mass", 100, 0., 8.);
  hIMAll->Sumw2();

  TH1F *hIMOppositeCharge =
      new TH1F("H", "Inv. Mass - opposite charge", 100, 0., 8.);
  hIMOppositeCharge->Sumw2();

  TH1F *hIMSameCharge = new TH1F("I", "Inv. Mass - same charge", 100, 0., 8.);
  hIMSameCharge->Sumw2();

  TH1F *hIMOppositePK =
      new TH1F("J", "Inv. Mass - opposite charge pions and kaons", 100, 0., 8.);
  hIMOppositePK->Sumw2();

  TH1F *hIMSamePK =
      new TH1F("K", "Inv. Mass - same charge pions and kaons", 100, 0., 8.);
  hIMSamePK->Sumw2();

  TH1F *hIMDecayParticles =
      new TH1F("L", "Inv. Mass - decay particles", 100, 0.7, 1.1);
  hIMDecayParticles->Sumw2();

  // inizio dei 1e5 eventi
  for (int n = 0; n < 1e5; ++n) {
    for (int i = 0; i < 100; ++i) {
      double theta = gRandom->Uniform(0, M_PI);
      hPolar->Fill(theta);

      double phi = gRandom->Uniform(0, 2 * M_PI);
      hAzimuthal->Fill(phi);

      double p = gRandom->Exp(1.);
      eventParticles[i].SetP(p * std::sin(theta) * std::cos(phi),
                             p * std::sin(theta) * std::sin(phi),
                             p * std::cos(theta));

      int random = gRandom->Uniform(1., 101.);
      if (random <= 80) {
        if (random <= 40) {
          eventParticles[i].SetIndex("pione+");
        } else {
          eventParticles[i].SetIndex("pione-");
        }
      } else if (random <= 90) {
        if (random <= 85) {
          eventParticles[i].SetIndex("Kaone+");
        } else {
          eventParticles[i].SetIndex("Kaone-");
        }
      } else if (random <= 99) {
        if (std::rand() % 2) {
          eventParticles[i].SetIndex("protone+");
        } else {
          eventParticles[i].SetIndex("protone-");
        }
      } else {
        eventParticles[i].SetIndex("K*");

        Particle p1{};
        Particle p2{};

        if (std::rand() % 2) {
          p1.SetIndex("pione+");
          p2.SetIndex("Kaone-");
        } else {
          p1.SetIndex("pione-");
          p2.SetIndex("Kaone+");
        }

        eventParticles[i].Decay2body(p1, p2);

        decayParticles.push_back(p1);
        decayParticles.push_back(p2);
      }
    }  // fine della generazione delle ~100 particelle

    eventParticles.insert(eventParticles.end(), decayParticles.begin(),
                          decayParticles.end());

    // riempimento istogrammi
    for (unsigned long i = 0; i < eventParticles.size(); ++i) {
      hParticleTypes->Fill(eventParticles[i].GetIndex());
      hImpulse->Fill(eventParticles[i].GetP());
      hTrsImpulse->Fill(eventParticles[i].GetTrsP());
      hEnergy->Fill(eventParticles[i].TotalEnergy());

      // istogrammi massa invariante
      for (unsigned long j = i + 1; j < eventParticles.size(); ++j) {
        // tutte le particelle
        hIMAll->Fill(eventParticles[i].InvariantMass(eventParticles[j]));

        // segno discorde
        if (eventParticles[i].GetCharge() * eventParticles[j].GetCharge() ==
            -1) {
          hIMOppositeCharge->Fill(
              eventParticles[i].InvariantMass(eventParticles[j]));

          // pioni e kaoni discordi
          if ((eventParticles[i].GetMass() + eventParticles[j].GetMass()) ==
              0.63324) {
            hIMOppositePK->Fill(
                eventParticles[i].InvariantMass(eventParticles[j]));
          }
        }
        // segno concorde
        else if (eventParticles[i].GetCharge() *
                     eventParticles[j].GetCharge() ==
                 1) {
          hIMSameCharge->Fill(
              eventParticles[i].InvariantMass(eventParticles[j]));

          // pioni e kaoni concordi
          if ((eventParticles[i].GetMass() + eventParticles[j].GetMass()) ==
              0.63324) {
            hIMSamePK->Fill(eventParticles[i].InvariantMass(eventParticles[j]));
          }
        }
      }
    }

    // particelle di decadimento
    for (unsigned long i = 0; i < decayParticles.size(); i += 2) {
      hIMDecayParticles->Fill(
          decayParticles[i].InvariantMass(decayParticles[i + 1]));
    }

    eventParticles.erase(eventParticles.end() - decayParticles.size(),
                         eventParticles.end());
    decayParticles.clear();

  }  // fine dei 1e5 eventi

  // inserimento degli istogrammi in un file ROOT
  TFile *file = new TFile("histograms.root", "RECREATE");

  hParticleTypes->Write();
  hPolar->Write();
  hAzimuthal->Write();
  hImpulse->Write();
  hTrsImpulse->Write();
  hEnergy->Write();
  hIMAll->Write();
  hIMOppositeCharge->Write();
  hIMSameCharge->Write();
  hIMOppositePK->Write();
  hIMSamePK->Write();
  hIMDecayParticles->Write();

  file->Close();

  hImpulse->Draw();

  return 0;
}
