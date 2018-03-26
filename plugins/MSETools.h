#ifndef MSETOOLS_H
#define MSETOOLS_H

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"




namespace mse
{
enum MuonParentage {PROMPT, HF, LF, NOT_A_MUON, NMU_PAR_TYPES};
extern const char* enumNames[4];

bool isBHadron (int pdgId);
bool decaysToB (const reco::GenParticle& gp);
bool isCHadron (int pdgId);
bool decaysToC (const reco::GenParticle& gp);

const pat::PackedGenParticle getMatchedGenParticle(const pat::GenericParticle &, const edm::View<pat::PackedGenParticle>&, int absPdgId);
//const pat::PackedGenParticle getMatchedGenParticle(const pat::GenericParticle &, const edm::View<pat::PackedGenParticle>&, int absPdgId);

const reco::GenParticle getMotherPacked(const pat::PackedGenParticle&);
const reco::GenParticle getMotherPacked(const reco::GenParticle&);

const reco::GenParticle getMother(const reco::GenParticle&);
//const reco::GenJet match(pat::Jet jet, double resolution,const reco::GenJetCollection m_genJets);

MuonParentage getParentType(const reco::GenParticle&);

bool isGoodVertex(const reco::Vertex&);
}

#endif
