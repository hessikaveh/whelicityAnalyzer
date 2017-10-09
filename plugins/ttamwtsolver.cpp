#include "ttamwtsolver.h"


void TtAMWTSolver::SetConstraints(const double xx, const double yy)
{
    pxmiss_ = xx;
    pymiss_ = yy;
}
string TtAMWTSolver::getEnvVar( string const & key ) const
{
    char * val = getenv( key.c_str() );
    return val == NULL ? string("") : string(val);
}
TtAMWTSolver::TtAMWTSolver(bool isData, const double b, const double e, const double s, const double mW, const double mB, const string ptRes, const string phiRes, const string sfRes):topmass_begin(b),
    topmass_end(e),
    topmass_step(s),
    mw(mW),
    mb(mB),
    pxmiss_(0),
    pymiss_(0),
    e_com(13000)

{
    if (isData) nbrJetSmear = 1000;
    else nbrJetSmear = 1000;
//    cout << "Number of iterations: " << nbrJetSmear<<endl;

    //PDF Initialization
    string path;
    path = getEnvVar("LHAPDF_DATA_PATH");
    string pdfSet = "NNPDF30_nlo_as_0118";
    LHAPDF::setPDFPath(path);
    LHAPDF::initPDFSet(pdfSet);

    //    string dataPath = getEnvVar("JER_DATA_PATH");
    //    string dataPath = "./data/";
    string ptFileName,  phiFileName,SFFileName;
    ptFileName = ptRes;
    phiFileName = phiRes;
    SFFileName = sfRes;

    //    ptFileName = "Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt";
    //    etaFileName = dataPath + "Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt";
    //    phiFileName = "Spring16_25nsV10_MC_PhiResolution_AK4PFchs.txt";
    //    SFFileName = dataPath + "Spring16_25nsV10_MC_SF_AK4PFchs.txt";
    //TRandom 3 class
    rand3 = new TRandom3();

    //    bool doGaussian = true;
//    cout <<ptFileName<<endl;
    ptResol  = new JME::JetResolution(ptFileName);

//    cout <<phiFileName<<endl;
    phiResol = new JME::JetResolution(phiFileName);

//    cout << SFFileName <<endl;
    SF = new JME::JetResolutionScaleFactor(SFFileName);
//    std::uint32_t seed = (std::uint32_t) 37428479;

    h_ptsm = new TH1F("h_ptsm","",100,0,300);
    h_etasm = new TH1F("h_etasm","",100,-3,3);
    h_phism = new TH1F("h_phism","",100,-5,5);
    h_SmearF = new TH1F("h_SmearFactorPt","",600,-3,3);
    h_SmearFPhi =new TH1F("h_SmearFactorPhi","",600,-3,3);
    h_SmearFEta =new TH1F("h_SmearFactorEta","",600,-3,3);
    h_Met = new TH2F("h_MET","",400,-400,400,400,-400,400);
    h_topMass = new TH1F("topMass","",100,172.4,173.1);
    f_out = TFile::Open("sanityCheckHistos.root","RECREATE");


}

TtAMWTSolver::NeutrinoSolution TtAMWTSolver::NuSolver(const TLorentzVector &LV_l, const TLorentzVector &LV_l_, const TLorentzVector &LV_b, const TLorentzVector &LV_b_)
{
    math::XYZTLorentzVector maxLV_n  = math::XYZTLorentzVector(0,0,0,0);
    math::XYZTLorentzVector maxLV_n_ = math::XYZTLorentzVector(0,0,0,0);

    //Shift the jets by an energy scale factor
    TVector2 met_sh;

    met_sh.Set(pxmiss_,pymiss_);
    vector<TLorentzVector> jets_by_pt_sh;
    for (int n_bjet_dummy =0; n_bjet_dummy < 2; ++n_bjet_dummy)
    {
//        cout << "pre " <<n_bjet_dummy << endl;

        if(n_bjet_dummy == 0)
        {
            jets_by_pt_sh.push_back(LV_b);
            jets_by_pt_sh.push_back(LV_b_);
//            cout << "1 " <<n_bjet_dummy << endl;
        }
        else if (n_bjet_dummy == 1)
        {
            jets_by_pt_sh.push_back(LV_b_);
            jets_by_pt_sh.push_back(LV_b);
//            cout << "2 " <<n_bjet_dummy << endl;

        }





        TVector2 met_sm;
        vector <TLorentzVector> jets_by_pt_sm;
        for (int iSM_JMT = 0; iSM_JMT < nbrJetSmear; iSM_JMT++)
        {


            if (nbrJetSmear == 1) dont_smear_JetMET(jets_by_pt_sh, met_sh, jets_by_pt_sm, met_sm);
            else smear_JetMET(jets_by_pt_sh, met_sh, jets_by_pt_sm, met_sm, rand3, (LV_l + LV_l_));
//        jet_energy_scale_factor.push_back(jets_by_pt_sm.at(0).Pt());
//        cout << jets_by_pt_sm.at(0).Pt() << endl;
            h_etasm->Fill(jets_by_pt_sm.at(0).Eta());
            h_ptsm->Fill(jets_by_pt_sm.at(0).Pt());
            h_phism->Fill(jets_by_pt_sm.at(0).Phi());
            h_SmearF->Fill(jets_by_pt_sm.at(0).Pt()/jets_by_pt_sh.at(0).Pt());
            h_SmearFPhi->Fill(jets_by_pt_sm.at(0).Phi()/jets_by_pt_sh.at(0).Phi());
            h_SmearFEta->Fill(jets_by_pt_sm.at(0).Eta()/jets_by_pt_sh.at(0).Eta());
            h_Met->Fill(met_sh.X(),met_sm.X());
            //loop on top mass parameter
            double weightmaxdum = -1;
            double mt = 172.5;
//            for(double mt = topmass_begin;
//                    mt < topmass_end + 0.5*topmass_step;
//                    mt += topmass_step)
            //  {
            double q_coeff[5], q_sol[4];
            // FindCoeff(LV_l, LV_l_, LV_b, LV_b_, mt, mt, pxmiss_, pymiss_, q_coeff);
            FindCoeff(LV_l, LV_l_, jets_by_pt_sm[0], jets_by_pt_sm[1], mt, mt, met_sm.X(), met_sm.Y(), q_coeff);
            int NSol = quartic(q_coeff, q_sol);

            //loop on all solutions
            for (int isol = 0; isol < NSol; isol++)
            {
                //TopRec(LV_l, LV_l_, LV_b, LV_b_, q_sol[isol]);
                TopRec(LV_l, LV_l_, jets_by_pt_sm[0], jets_by_pt_sm[1], q_sol[isol]);
                double weight = get_weight(jets_by_pt_sm[0],jets_by_pt_sm[1],LV_l,LV_l_,LV_n,LV_n_,mt);
//		cout <<"weight of mt answer: " <<  weight << endl;
                h_topMass->Fill(mt,weight);
                if (weight > weightmaxdum)
                {
                    weightmaxdum =weight;
                    maxLV_n.SetPxPyPzE(LV_n.Px(), LV_n.Py(), LV_n.Pz(), LV_n.E());
                    maxLV_n_.SetPxPyPzE(LV_n_.Px(), LV_n_.Py(), LV_n_.Pz(), LV_n_.E());
                    weightmax = weightmaxdum;
                }
            }
            //}
        }//End of JetMet Smearing loop
    }
    TtAMWTSolver::NeutrinoSolution nuSol;
    nuSol.neutrino    = reco::LeafCandidate(0, maxLV_n  );
    nuSol.neutrinoBar = reco::LeafCandidate(0, maxLV_n_ );
    nuSol.weight = weightmax;



    return nuSol;
}






void TtAMWTSolver::writeOut()
{
    f_out->cd();
    h_ptsm->Write();
    h_phism->Write();
    h_SmearF->Write();
    h_etasm->Write();
    h_Met->Write();
    h_SmearFEta->Write();
    h_SmearFPhi->Write();
    h_topMass->Write();
    f_out->Close();
}

double TtAMWTSolver::get_weight(const TLorentzVector & bquark1, const TLorentzVector & bquark2,
                                const TLorentzVector & lep_p, const TLorentzVector & lep_m, const TLorentzVector & nu1,
                                const TLorentzVector & nu2, const double top_mass) const
{

    TLorentzVector t1 = lep_p + nu1 + bquark1;
    TLorentzVector t2 = lep_m + nu2 + bquark2;

    //Get
    double prob_dalitz = 1.0;
    prob_dalitz *= get_dalitz_prob( lep_p, t1, mb, mw );
    prob_dalitz *= get_dalitz_prob( lep_m, t2, mb, mw );

    //Determine x1 and x2
    double x1 = ( t1.E() + t2.E() + t1.Pz() + t2.Pz() ) / e_com;
    double x2 = ( t1.E() + t2.E() - t1.Pz() - t2.Pz() ) / e_com;

    vector <double>  f1 = LHAPDF::xfx(x1, top_mass);
    vector <double>  f2 = LHAPDF::xfx(x2, top_mass);

    // The order of f:
    //    -t  -b  -c  -s  -u  -d   g   d   u   s   c   b   t
    //    -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6
    //     0   1   2   3   4   5   6   7   8   9   10  11  12

    double sbar1 = f1[3], sbar2 = f2[3];
    double ubar1 = f1[4], ubar2 = f2[4];
    double dbar1 = f1[5], dbar2 = f2[5];
    double g1    = f1[6], g2    = f2[6];
    double d1    = f1[7], d2    = f2[7];
    double u1    = f1[8], u2    = f2[8];
    double s1    = f1[9], s2    = f2[9];

    //Should glue-glue be doubled? Probably not, but plot histo later
    double pdf_prob = (u1*ubar2 + u2*ubar1 +
                       d1*dbar2 + d2*dbar1 +
                       s1*sbar2 + s2*sbar1 +
                       g1*g2);

    //   double tt_pt_prob = get_top_pt_prob(t1.Pt()) * get_top_pt_prob(t2.Pt());

    //   double two_bjet_prob = get_2bjet_prob(bquark1, bquark2, mapJetPhi2Discr);

    //   double s_weight = pdf_prob * prob_dalitz * tt_pt_prob * two_bjet_prob;
    //   double s_weight = pdf_prob * prob_dalitz * two_bjet_prob;
    //   double s_weight = pdf_prob * prob_dalitz * tt_pt_prob;
    return pdf_prob * prob_dalitz;
}

double TtAMWTSolver::get_dalitz_prob(const TLorentzVector & lep, const TLorentzVector & top,
                                     double mb, double mw) const
{
    double mte = lep.Dot( top );
    double mt = top.M();
    double mt2 = mt * mt;
    double mb2 = mb * mb;
    double mw2 = mw * mw;
    double mt2_mb2 = mt2 - mb2;

    return 4. * mte * ( mt2 - mb2 - 2. * mte ) /
           ( mt2_mb2 * mt2_mb2 + mw2 * ( mt2 + mb2 ) - 2. * mw2 * mw2 );
}

double TtAMWTSolver::get_2bjet_prob(const TLorentzVector & jet1, const TLorentzVector & jet2,
                                    map <double,double> & mapJetPhi2Discr) const
{

    double proba1 =exp(-4*mapJetPhi2Discr[jet1.Phi()]);
    double proba2 =exp(-4*mapJetPhi2Discr[jet2.Phi()]);

    return (1-proba1)*(1-proba2);
}


double TtAMWTSolver::get_top_pt_prob(const double pt) const
{

    double prob_array[120] =
    {
        0.00200245,	 0.0060801,	0.00999056,	0.013597,	0.0170136,	0.0203661,	0.0231052,	0.0259027,	0.0279369,	0.0299741,
        0.0313266,	 0.0328515,	0.0333733,	0.0339433,	0.0341794,	0.0340638,	0.0339737,	0.0336051,	0.0332196,	0.0317563,
        0.0308568,	 0.0299114,	0.028731,	0.0275211,	0.0260246,	0.0246938,	0.0234093,	0.022064,	0.0207067,	0.0193499,
        0.0179213,	 0.0169745,	0.0156669,	0.0145741,	0.013649,	0.012431,	0.0114861,	0.0106781,	0.00984938,	0.00909871,
        0.00843331,	 0.00785657,	0.00723358,	0.0066419,	0.00599867,	0.00563297,	0.00518391,	0.004678,	0.00437831,	0.00410608,
        0.00364064,	 0.00340744,	0.00303836,	0.00285093,	0.0026134,	0.00236285,	0.00221879,	0.00197884,	0.00180586,	0.00173118,
        0.00154375,	 0.00147389,	0.00136018,	0.0012378,	0.00110433,	0.0010099,	0.000963642,	0.000868241,	0.000836923,	0.000752604,
        0.000649495, 0.000618176,	0.000603722,	0.000543976,	0.000476521,	0.000447612,	0.000396057,	0.000378229,	0.000363293,	0.000354138,
        0.000287647, 0.000282347,	0.00025922,	0.000260183,	0.000238983,	0.000201883,	0.000204292,	0.000198992,	0.000176828,	0.000161892,
        0.00013491,	 0.00014551,	0.000106964,	0.000112264,	0.000101664,	0.00010311,	9.29914e-05,	7.8055e-05,	7.8055e-05,	6.69731e-05,
        7.08277e-05, 5.49276e-05,	6.07094e-05,	6.79367e-05,	4.77003e-05,	4.77003e-05,	4.19184e-05,	3.71002e-05,	4.43275e-05,	3.13184e-05,
        2.98729e-05, 3.66184e-05,	2.60183e-05,	2.89093e-05,	2.55365e-05,	2.55365e-05,	2.45729e-05,	2.12001e-05,	2.16819e-05,	2.65001e-05
    };

    //The array is binned in intervals of 5 GeV
    int ipt = (int) pt / 5;

    double pt_prob = 0;

    if (ipt > 119) pt_prob = 2E-5;
    else           pt_prob = prob_array[ipt];

    return pt_prob;
}




void TtAMWTSolver::smear_JetMET(const vector <TLorentzVector> & orig_jets, const TVector2 & orig_met,
                                vector <TLorentzVector> & smear_jets, TVector2 & smear_met,
                                TRandom3* rand3,  const TLorentzVector & lep_sum) const
{

    smear_jets.clear();

    double sum_jpx = 0;
    double sum_jpy = 0;

    double sum_jpx_sm = 0;
    double sum_jpy_sm = 0;

    double Pt_sm, Eta_sm, Phi_sm;

    TLorentzVector v_temp;
    std::mt19937_64 m_random_generator(21213232546);
    // Seed with a real random value, if available
    std::random_device r;

    // Choose a random mean between 1 and 6
//        std::default_random_engine e1(r());
//        std::uniform_int_distribution<int> uniform_dist(1, 6);
//        int mean = uniform_dist(e1);
//        std::cout << "Randomly-chosen mean: " << mean << '\n';

    // Generate a normal distribution around that mean
    std::seed_seq seed2{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 e2(seed2);
//        std::normal_distribution<> normal_dist(mean, 2);



    double smearFactorPt = 1.;
    double smearFactorPhi = 1.;
    for (unsigned int sui = 0; sui < orig_jets.size(); sui++)
    {

        JME::JetParameters par;
        par.setJetEta(orig_jets.at(sui).Eta());
        par.setJetPt(orig_jets.at(sui).Pt());
        par.setRho(rho);
        double jet_resPt = ptResol->getResolution(par);
        double jet_resPhi = phiResol->getResolution(par);
        double jer_sf = SF->getScaleFactor({{JME::Binning::JetEta, orig_jets.at(sui).Eta()}},Variation::NOMINAL);
        double sigma = jet_resPt;// * std::sqrt(jer_sf * jer_sf - 1);
        double sigma2 = jet_resPhi;
//        cout<< sigma <<endl;
        std::normal_distribution<> d(0, sigma);
        std::normal_distribution<> l(0,sigma2);
//        cout << d(e2)<< endl;
        smearFactorPt = 1. + d(e2)*std::sqrt(std::max(jer_sf*jer_sf -1,0.));
        smearFactorPhi = 1.+ l(e2);

        Pt_sm = orig_jets.at(sui).Pt()*smearFactorPt;
        Eta_sm = orig_jets.at(sui).Eta();
        Phi_sm = orig_jets.at(sui).Phi()*smearFactorPhi;
        //        cout <<"Jet "<< smearFactor<<" " << Pt_sm<<" " << Eta_sm<<" " << Phi_sm << endl;
        v_temp.SetPtEtaPhiE(Pt_sm, Eta_sm, Phi_sm, orig_jets.at(sui).E());
        //v_temp = orig_jets[sui]*smearFactorPt;
        sum_jpx += orig_jets.at(sui).Px();
        sum_jpy += orig_jets.at(sui).Py();

        sum_jpx_sm += v_temp.Px();
        sum_jpy_sm += v_temp.Py();

        smear_jets.push_back(v_temp);
//        jet_energy_scale_factor.push_back(smearFactor);
    }

    double unclust_metx = orig_met.Px() + sum_jpx + lep_sum.Px();
    double unclust_mety = orig_met.Py() + sum_jpy + lep_sum.Py();
//    h_SmearF->Fill(jet_energy_scale_factor[0]);
//    h_SmearF->Fill(jet_energy_scale_factor[1]);

    //10% resolution
    std::normal_distribution<> t(0, 0.10);
    double MetSmear = t(e2);
    double unclust_metx_sm = unclust_metx * (1 + MetSmear );
    double unclust_mety_sm = unclust_mety * (1 + MetSmear );
//    cout << unclust_metx_sm<<"Met " << unclust_mety_sm << endl;

    smear_met.Set(orig_met.Px() + sum_jpx - unclust_metx - sum_jpx_sm + unclust_metx_sm, orig_met.Py() + sum_jpy - unclust_mety - sum_jpy_sm + unclust_mety_sm);
//smear_met.Set(unclust_metx_sm -sum_jpx + sum_jpx_sm, unclust_mety_sm - sum_jpy + sum_jpy_sm);
}

void TtAMWTSolver::dont_smear_JetMET(const vector <TLorentzVector> & orig_jets,
                                     const TVector2 & orig_met, vector <TLorentzVector> & smear_jets,
                                     TVector2 & smear_met) const
{
    smear_met  = orig_met;
    smear_jets = orig_jets;
}

void TtAMWTSolver::FindCoeff(const TLorentzVector &al, const TLorentzVector &l, const TLorentzVector &b_al, const TLorentzVector &b_l, const double mt, const double mat, const double px_miss, const double py_miss, double* koeficienty)
{
    double E, apom1, apom2, apom3;
    double k11, k21, k31, k41,  cpom1, cpom2, cpom3, l11, l21, l31, l41, l51, l61, k1, k2, k3, k4, k5,k6;
    double l1, l2, l3, l4, l5, l6, k15, k25, k35, k45;

    C = -al.Px()-b_al.Px()-l.Px()-b_l.Px() + px_miss;
    D = -al.Py()-b_al.Py()-l.Py()-b_l.Py() + py_miss;

    // right side of first two linear equations - missing pT

    E = (sqr(mt)-sqr(mw)-sqr(mb))/(2*b_al.E())-sqr(mw)/(2*al.E())-al.E()+al.Px()*b_al.Px()/b_al.E()+al.Py()*b_al.Py()/b_al.E()+al.Pz()*b_al.Pz()/b_al.E();
    F = (sqr(mat)-sqr(mw)-sqr(mb))/(2*b_l.E())-sqr(mw)/(2*l.E())-l.E()+l.Px()*b_l.Px()/b_l.E()+l.Py()*b_l.Py()/b_l.E()+l.Pz()*b_l.Pz()/b_l.E();

    m1 = al.Px()/al.E()-b_al.Px()/b_al.E();
    m2 = al.Py()/al.E()-b_al.Py()/b_al.E();
    m3 = al.Pz()/al.E()-b_al.Pz()/b_al.E();

    n1 = l.Px()/l.E()-b_l.Px()/b_l.E();
    n2 = l.Py()/l.E()-b_l.Py()/b_l.E();
    n3 = l.Pz()/l.E()-b_l.Pz()/b_l.E();

    pom = E-m1*C-m2*D;
    apom1 = sqr(al.Px())-sqr(al.E());
    apom2 = sqr(al.Py())-sqr(al.E());
    apom3 = sqr(al.Pz())-sqr(al.E());

    k11 = 1/sqr(al.E())*(pow(mw,4)/4+sqr(C)*apom1+sqr(D)*apom2+apom3*sqr(pom)/sqr(m3)+sqr(mw)*(al.Px()*C+al.Py()*D+al.Pz()*pom/m3)+2*al.Px()*al.Py()*C*D+2*al.Px()*al.Pz()*C*pom/m3+2*al.Py()*al.Pz()*D*pom/m3);
    k21 = 1/sqr(al.E())*(-2*C*m3*n3*apom1+2*apom3*n3*m1*pom/m3-sqr(mw)*m3*n3*al.Px()+sqr(mw)*m1*n3*al.Pz()-2*al.Px()*al.Py()*D*m3*n3+2*al.Px()*al.Pz()*C*m1*n3-2*al.Px()*al.Pz()*n3*pom+2*al.Py()*al.Pz()*D*m1*n3);
    k31 = 1/sqr(al.E())*(-2*D*m3*n3*apom2+2*apom3*n3*m2*pom/m3-sqr(mw)*m3*n3*al.Py()+sqr(mw)*m2*n3*al.Pz()-2*al.Px()*al.Py()*C*m3*n3+2*al.Px()*al.Pz()*C*m2*n3-2*al.Py()*al.Pz()*n3*pom+2*al.Py()*al.Pz()*D*m2*n3);
    k41 = 1/sqr(al.E())*(2*apom3*m1*m2*sqr(n3)+2*al.Px()*al.Py()*sqr(m3)*sqr(n3)-2*al.Px()*al.Pz()*m2*m3*sqr(n3)-2*al.Py()*al.Pz()*m1*m3*sqr(n3));
    k51 = 1/sqr(al.E())*(apom1*sqr(m3)*sqr(n3)+apom3*sqr(m1)*sqr(n3)-2*al.Px()*al.Pz()*m1*m3*sqr(n3));
    k61 = 1/sqr(al.E())*(apom2*sqr(m3)*sqr(n3)+apom3*sqr(m2)*sqr(n3)-2*al.Py()*al.Pz()*m2*m3*sqr(n3));

    cpom1 = sqr(l.Px())-sqr(l.E());
    cpom2 = sqr(l.Py())-sqr(l.E());
    cpom3 = sqr(l.Pz())-sqr(l.E());

    l11 = 1/sqr(l.E())*(pow(mw,4)/4+cpom3*sqr(F)/sqr(n3)+sqr(mw)*l.Pz()*F/n3);
    l21 = 1/sqr(l.E())*(-2*cpom3*F*m3*n1/n3+sqr(mw)*(l.Px()*m3*n3-l.Pz()*n1*m3)+2*l.Px()*l.Pz()*F*m3);
    l31 = 1/sqr(l.E())*(-2*cpom3*F*m3*n2/n3+sqr(mw)*(l.Py()*m3*n3-l.Pz()*n2*m3)+2*l.Py()*l.Pz()*F*m3);
    l41 = 1/sqr(l.E())*(2*cpom3*n1*n2*sqr(m3)+2*l.Px()*l.Py()*sqr(m3)*sqr(n3)-2*l.Px()*l.Pz()*n2*n3*sqr(m3)-2*l.Py()*l.Pz()*n1*n3*sqr(m3));
    l51 = 1/sqr(l.E())*(cpom1*sqr(m3)*sqr(n3)+cpom3*sqr(n1)*sqr(m3)-2*l.Px()*l.Pz()*n1*n3*sqr(m3));
    l61 = 1/sqr(l.E())*(cpom2*sqr(m3)*sqr(n3)+cpom3*sqr(n2)*sqr(m3)-2*l.Py()*l.Pz()*n2*n3*sqr(m3));

    k1 = k11*k61;
    k2 = k61*k21/k51;
    k3 = k31;
    k4 = k41/k51;
    k5 = k61/k51;
    k6 = 1;

    l1 = l11*k61;
    l2 = l21*k61/k51;
    l3 = l31;
    l4 = l41/k51;
    l5 = l51*k61/(sqr(k51));
    l6 = l61/k61;

    k15 = k1*l5-l1*k5;
    k25 = k2*l5-l2*k5;
    k35 = k3*l5-l3*k5;
    k45 = k4*l5-l4*k5;

    k16 = k1*l6-l1*k6;
    k26 = k2*l6-l2*k6;
    k36 = k3*l6-l3*k6;
    k46 = k4*l6-l4*k6;
    k56 = k5*l6-l5*k6;

    koeficienty[0] = k15*sqr(k36)-k35*k36*k16-k56*sqr(k16);
    koeficienty[1] = 2*k15*k36*k46+k25*sqr(k36)+k35*(-k46*k16-k36*k26)-k45*k36*k16-2*k56*k26*k16;
    koeficienty[2] = k15*sqr(k46)+2*k25*k36*k46+k35*(-k46*k26-k36*k56)-k56*(sqr(k26)+2*k56*k16)-k45*(k46*k16+k36*k26);
    koeficienty[3] = k25*sqr(k46)-k35*k46*k56-k45*(k46*k26+k36*k56)-2*sqr(k56)*k26;
    koeficienty[4] = -k45*k46*k56-pow(k56,3);

    // normalization of coefficients
    int moc=(int(log10(fabs(koeficienty[0])))+int(log10(fabs(koeficienty[4]))))/2;

    koeficienty[0]=koeficienty[0]/TMath::Power(10,moc);
    koeficienty[1]=koeficienty[1]/TMath::Power(10,moc);
    koeficienty[2]=koeficienty[2]/TMath::Power(10,moc);
    koeficienty[3]=koeficienty[3]/TMath::Power(10,moc);
    koeficienty[4]=koeficienty[4]/TMath::Power(10,moc);
}

void TtAMWTSolver::TopRec(const TLorentzVector &al, const TLorentzVector &l, const TLorentzVector &b_al, const TLorentzVector &b_l, const double sol)
{
    TVector3 t_ttboost;
    TLorentzVector aux;
    double pxp, pyp, pzp, pup, pvp, pwp;

    pxp = sol*(m3*n3/k51);
    pyp = -(m3*n3/k61)*(k56*pow(sol,2) + k26*sol + k16)/(k36 + k46*sol);
    pzp = -1/n3*(n1*pxp + n2*pyp - F);
    pwp = 1/m3*(m1*pxp + m2*pyp + pom);
    pup = C - pxp;
    pvp = D - pyp;

    LV_n_.SetXYZM(pxp, pyp, pzp, 0.0);
    LV_n.SetXYZM(pup, pvp, pwp, 0.0);

    LV_t_ = b_l + l + LV_n_;
    LV_t = b_al + al + LV_n;

    aux = (LV_t_ + LV_t);
    t_ttboost = -aux.BoostVector();
    LV_tt_t_ = LV_t_;
    LV_tt_t = LV_t;
    LV_tt_t_.Boost(t_ttboost);
    LV_tt_t.Boost(t_ttboost);
}

int TtAMWTSolver::quartic(double *koeficienty, double *koreny) const
{
    double w, b0, b1, b2;
    double c[4];
    double d0, d1, h, t, z;
    double *px;

    if (koeficienty[4]==0.0)
        return cubic(koeficienty, koreny);
    /* quartic problem? */
    w = koeficienty[3]/(4*koeficienty[4]);
    /* offset */
    b2 = -6*sqr(w) + koeficienty[2]/koeficienty[4];
    /* koeficienty. of shifted polynomial */
    b1 = (8*sqr(w) - 2*koeficienty[2]/koeficienty[4])*w + koeficienty[1]/koeficienty[4];
    b0 = ((-3*sqr(w) + koeficienty[2]/koeficienty[4])*w - koeficienty[1]/koeficienty[4])*w + koeficienty[0]/koeficienty[4];

    c[3] = 1.0;
    /* cubic resolvent */
    c[2] = b2;
    c[1] = -4*b0;
    c[0] = sqr(b1) - 4*b0*b2;

    cubic(c, koreny);
    z = koreny[0];
    //double z1=1.0,z2=2.0,z3=3.0;
    //TMath::RootsCubic(c,z1,z2,z3);
    //if (z2 !=0) z = z2;
    //if (z1 !=0) z = z1;
    /* only lowermost root needed */

    int nreal = 0;
    px = koreny;
    t = sqrt(0.25*sqr(z) - b0);
    for(int i=-1; i<=1; i+=2)
    {
        d0 = -0.5*z + i*t;
        /* coeffs. of quadratic factor */
        d1 = (t!=0.0)? -i*0.5*b1/t : i*sqrt(-z - b2);
        h = 0.25*sqr(d1) - d0;
        if (h>=0.0)
        {
            h = sqrt(h);
            nreal += 2;
            *px++ = -0.5*d1 - h - w;
            *px++ = -0.5*d1 + h - w;
        }
    }

    //  if (nreal==4) {
    /* sort results */
    //    if (koreny[2]<koreny[0]) SWAP(koreny[0], koreny[2]);
    //    if (koreny[3]<koreny[1]) SWAP(koreny[1], koreny[3]);
    //    if (koreny[1]<koreny[0]) SWAP(koreny[0], koreny[1]);
    //    if (koreny[3]<koreny[2]) SWAP(koreny[2], koreny[3]);
    //    if (koreny[2]<koreny[1]) SWAP(koreny[1], koreny[2]);
    //  }
    return nreal;
}

int TtAMWTSolver::cubic(const double *coeffs, double *koreny) const
{
    unsigned nreal;
    double w, p, q, dis, h, phi;

    if (coeffs[3]!=0.0)
    {
        /* cubic problem? */
        w = coeffs[2]/(3*coeffs[3]);
        p = sqr(coeffs[1]/(3*coeffs[3])-sqr(w))*(coeffs[1]/(3*coeffs[3])-sqr(w));
        q = -0.5*(2*sqr(w)*w-(coeffs[1]*w-coeffs[0])/coeffs[3]);
        dis = sqr(q)+p;
        /* discriminant */
        if (dis<0.0)
        {
            /* 3 real solutions */
            h = q/sqrt(-p);
            if (h>1.0) h = 1.0;
            /* confine the argument of */
            if (h<-1.0) h = -1.0;
            /* acos to [-1;+1] */
            phi = acos(h);
            p = 2*TMath::Power(-p, 1.0/6.0);
            for(unsigned i=0; i<3; i++)
                koreny[i] = p*cos((phi+2*i*TMath::Pi())/3.0) - w;
            if (koreny[1]<koreny[0]) SWAP(koreny[0], koreny[1]);
            /* sort results */
            if (koreny[2]<koreny[1]) SWAP(koreny[1], koreny[2]);
            if (koreny[1]<koreny[0]) SWAP(koreny[0], koreny[1]);
            nreal = 3;
        }
        else
        {
            /* only one real solution */
            dis = sqrt(dis);
            h = TMath::Power(fabs(q+dis), 1.0/3.0);
            p = TMath::Power(fabs(q-dis), 1.0/3.0);
            koreny[0] = ((q+dis>0.0)? h : -h) + ((q-dis>0.0)? p : -p) -  w;
            nreal = 1;
        }

        /* Perform one step of a Newton iteration in order to minimize
         round-off errors */
        for(unsigned i=0; i<nreal; i++)
        {
            h = coeffs[1] + koreny[i] * (2 * coeffs[2] + 3 * koreny[i] * coeffs[3]);
            if (h != 0.0)
                koreny[i] -= (coeffs[0] + koreny[i] * (coeffs[1] + koreny[i] * (coeffs[2] + koreny[i] * coeffs[3])))/h;
        }
    }

    else if (coeffs[2]!=0.0)
    {
        /* quadratic problem? */
        p = 0.5*coeffs[1]/coeffs[2];
        dis = sqr(p) - coeffs[0]/coeffs[2];
        if (dis>=0.0)
        {
            /* two real solutions */
            dis = sqrt(dis);
            koreny[0] = -p - dis;
            koreny[1] = -p + dis;
            nreal = 2;
        }
        else
            /* no real solution */
            nreal = 0;
    }

    else if (coeffs[1]!=0.0)
    {
        /* linear problem? */
        koreny[0] = -coeffs[0]/coeffs[1];
        nreal = 1;
    }

    else
        /* no equation */
        nreal = 0;

    return nreal;
}

void TtAMWTSolver::SWAP(double &realone, double &realtwo) const
{
    if (realtwo < realone)
    {
        double aux = realtwo;
        realtwo = realone;
        realone = aux;
    }
}

TtAMWTSolver::~TtAMWTSolver()
{

}
