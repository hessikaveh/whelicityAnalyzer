#ifndef TTAMWTSOLVER_H
#define TTAMWTSOLVER_H
#include "DataFormats/Candidate/interface/LeafCandidate.h"
//#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include <memory>
#include <vector>
#include <sstream>
#include <iostream>
#include <random>
#include <time.h>
#include <TRandom3.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include "TLorentzVector.h"
#include "TMath.h"
#include <LHAPDF/LHAPDF.h>
#include "NeutrinoEllipseCalculator.h"

using namespace std;

class TtAMWTSolver
{

public:

    struct NeutrinoSolution {
        double weight;
        reco::LeafCandidate neutrino;
        reco::LeafCandidate neutrinoBar;
    };

    void SetConstraints(const double xx=0, const double yy=0);
    void SetRho(const double rho_){rho=rho_;}
    TtAMWTSolver(bool isData, const double b, const double e, const double s, const double mW, const double mB, const string ptRes,const string phiRes,const string sfRes);
    ~TtAMWTSolver();
    NeutrinoSolution NuSolver(const TLorentzVector &LV_l, const TLorentzVector &LV_l_, const TLorentzVector &LV_b, const TLorentzVector &LV_b_);



    void writeOut();

private:

    int nbrJetSmear;
    double get_weight(const TLorentzVector &bquark1, const TLorentzVector &bquark2, const TLorentzVector &lep_p, const TLorentzVector &lep_m, const TLorentzVector &nu1, const TLorentzVector &nu2, const double top_mass) const;
    double get_dalitz_prob(const TLorentzVector &lep, const TLorentzVector &top, double mb, double mw) const;
    void smear_JetMET(const vector<TLorentzVector> &orig_jets, const TVector2 &orig_met, vector<TLorentzVector> &smear_jets, TVector2 &smear_met, TRandom3 *rand3, const TLorentzVector &lep_sum) const;
    void dont_smear_JetMET(const vector<TLorentzVector> &orig_jets, const TVector2 &orig_met, vector<TLorentzVector> &smear_jets, TVector2 &smear_met) const;
    ///
    void FindCoeff(const TLorentzVector& al,
                   const TLorentzVector& l,
                   const TLorentzVector& b_al,
                   const TLorentzVector& b_l,
                   const double mt, const double mat, const double px_miss, const double py_miss,
                   double* koeficienty);
    ///
    void TopRec(const TLorentzVector& al,
                const TLorentzVector& l,
                const TLorentzVector& b_al,
                const TLorentzVector& b_l, const double sol);

    ///
    int quartic(double* koeficienty, double* koreny) const;
    ///
    int cubic(const double* coeffs, double* koreny) const;
    ///
    double sqr(const double x) const {return (x*x);}
    ///
    void SWAP(double& realone, double& realtwo) const;

private:
    ///
    const double topmass_begin;
    ///
    const double topmass_end;
    ///
    const double topmass_step;
    ///
     double mw;
    ///
    const double mb;
    ///
    double pxmiss_, pymiss_;
    const double e_com;

    double C;
    double D;
    double F;
    double pom;
    double k16;
    double k26;
    double k36;
    double k46;
    double k56;
    double k51;
    double k61;
    double m1;
    double m2;
    double m3;
    double n1;
    double n2;
    double n3;
    TH1F *h_ptsm;
    TH1F *h_etasm;
    TH1F *h_phism;
    TH1F *h_SmearF;
    TH1F *h_SmearFPhi;
    TH1F *h_SmearFEta;
    TH2F *h_Met;
    TH1F *h_topMass;
    ///
    TLorentzVector LV_n, LV_n_, LV_t, LV_t_, LV_tt_t, LV_tt_t_;
    /// provisional
    TLorentzVector genLV_n, genLV_n_;
//    std::mt19937 m_random_generator(std::uint32_t 87847938);
    TRandom3* rand3;
//    JetResolution* ptResol;
//    JetResolution* etaResol;
//    JetResolution* phiResol;
    double rho;
    JME::JetResolution* ptResol;

    JME::JetResolution* phiResol;
    JME::JetResolutionScaleFactor* SF;
    vector<double> jet_energy_scale_factor;

    double weightmax=0;

    double get_top_pt_prob(const double pt) const;
    double get_2bjet_prob(const TLorentzVector &jet1, const TLorentzVector &jet2, map<double, double> &mapJetPhi2Discr) const;
    string getEnvVar(const string &key) const;
    double Cofactor(TMatrixD M,int i, int j);
    TMatrixD factor_degenrate(TMatrixD G,int zero=0);
    TFile* f_out;
   // TH1F* h_WMASS;
  //  TH1F* h_WMASS2;

//    TFile* f_WMASS;
};

#endif // TTAMWTSOLVER_H
