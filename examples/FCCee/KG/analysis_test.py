#Mandatory: List of processes
processList = {
    # 1M b & c, 3M uds
    #'p8_ee_Zbb_ecm91':{'fraction':0.001, 'chunks':10}, #Run 0.1% statistics in two files named <outputDir>/p8_ee_Zbb_ecm91/chunk<N>.root
    #'p8_ee_Zcc_ecm91':{'fraction':0.001, 'chunks':10}, #Run 0.1% statistics in two files named <outputDir>/p8_ee_Zcc_ecm91/chunk<N>.root
    #'p8_ee_Zuds_ecm91':{'fraction':0.003, 'chunks':30}, #Run 0.3% statistics in one output file named <outputDir>/p8_ee_Zuds_ecm91.root

    # 400K b & c, 1M uds
    #'p8_ee_Zbb_ecm91':{'fraction':0.0004, 'chunks':4},
    #'p8_ee_Zcc_ecm91':{'fraction':0.0004, 'chunks':4},
    #'p8_ee_Zuds_ecm91':{'fraction':0.001, 'chunks':10},
    'p8_ee_Zuds_ecm91':{'fraction':0.0001},
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/spring2021/IDEA/"

#Optional: output directory, default is local running directory
#outputDir   = "outputs/FCCee/KG/"
outputDir   = "/afs/cern.ch/work/k/kgautam/private/latest/FCCAnalyses/outputs/FCCee/KG/"

#Optional: ncpus, default is 4
nCPUS       = 8

#Optional running on HTCondor, default is False
#runBatch    = True

#Optional batch queue name when running on HTCondor, default is workday
#batchQueue = "longlunch"

#Optional computing account when running on HTCondor, default is group_u_FCC.local_gen
#compGroup = "group_u_FCC.local_gen"

#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df2 = (
            df

            # alias
            .Alias("Particle0", "Particle#0.index") 
            .Alias("Particle1", "Particle#1.index") 
            .Alias("Jet3","Jet#3.index")
            .Alias("Muon0","Muon#0.index")
            .Alias("Electron0","Electron#0.index")
            .Alias("Photon0","Photon#0.index")
            .Define("Particle1_indices", "Particle1")

            # Reco'd particles
            .Define("RP_px",        "ReconstructedParticle::get_px(ReconstructedParticles)")
            .Define("RP_py",        "ReconstructedParticle::get_py(ReconstructedParticles)")
            .Define("RP_pz",        "ReconstructedParticle::get_pz(ReconstructedParticles)")               
            .Define("RP_m",         "ReconstructedParticle::get_mass(ReconstructedParticles)")
            .Define("RP_charge",    "ReconstructedParticle::get_charge(ReconstructedParticles)")
            #.Define("RP_hasTRK",    "ReconstructedParticle2Track::hasTRK(ReconstructedParticles)")
            .Define("RP_hasTRK",    "ROOT::VecOps::RVec<int> result; for(auto& charge : RP_charge){if(std::abs(charge)>0) result.push_back(1); else result.push_back(0);} return result;")
            .Define("RP_charged",   "ReconstructedParticle::sel_tag(1)(RP_hasTRK, ReconstructedParticles)")
            .Define("RP_neutral",   "ReconstructedParticle::sel_tag(0)(RP_hasTRK, ReconstructedParticles)")
            .Define("RP_isMuon",    "ReconstructedParticle::is_particle(Muon0, ReconstructedParticles)")
            .Define("RP_isElectron","ReconstructedParticle::is_particle(Electron0, ReconstructedParticles)")
            .Define("RP_isPhoton",  "ReconstructedParticle::is_particle(Photon0, ReconstructedParticles)")

            #Get charged PF vars 
            .Define("RP_charged_charge",    "ReconstructedParticle::get_charge(RP_charged)")
            .Define("RP_charged_p",         "ReconstructedParticle::get_p(RP_charged)")
            #.Define("RP_charged_p_template","ReconstructedParticle::get_p(RP_template)")
            .Define("RP_charged_theta",     "ReconstructedParticle::get_theta(RP_charged)")
            .Define("RP_charged_phi",       "ReconstructedParticle::get_phi(RP_charged)")
            .Define("RP_charged_mass",      "ReconstructedParticle::get_mass(RP_charged)") #should think about this, do we really mean mass of RP or assuming pion?
            .Define("RP_charged_Z0",        "ReconstructedParticle2Track::getRP2TRK_Z0(RP_charged, EFlowTrack_1)")
            .Define("RP_charged_D0",        "ReconstructedParticle2Track::getRP2TRK_D0(RP_charged, EFlowTrack_1)")
            .Define("RP_charged_Z0_sig",    "ReconstructedParticle2Track::getRP2TRK_Z0_sig(RP_charged, EFlowTrack_1)")
            .Define("RP_charged_D0_sig",    "ReconstructedParticle2Track::getRP2TRK_D0_sig(RP_charged, EFlowTrack_1)")

            #Get neutral PF vars 
            .Define("RP_neutral_p", "ReconstructedParticle::get_p(RP_neutral)")

            # First, reconstruct a vertex from all tracks
            .Define("VertexObject_allTracks",  "VertexFitterSimple::VertexFitter_Tk ( 1, EFlowTrack_1, true, 4.5, 20e-3, 300)")
            # Select the tracks that are reconstructed  as primaries
            .Define("RecoedPrimaryTracks",  "VertexFitterSimple::get_PrimaryTracks( VertexObject_allTracks, EFlowTrack_1, true, 4.5, 20e-3, 300, 0., 0., 0., 0)")
            .Define("n_RecoedPrimaryTracks",  "ReconstructedParticle2Track::getTK_n( RecoedPrimaryTracks )")
            # the final primary vertex : 
            .Define("PV", "VertexFitterSimple::VertexFitter_Tk ( 1, RecoedPrimaryTracks, true, 4.5, 20e-3, 300) ") # primary vertex object
            .Define("PrimaryVertex", "VertexingUtils::get_VertexData( PV )")
            
            # the secondary tracks
            .Define("SecondaryTracks",   "VertexFitterSimple::get_NonPrimaryTracks( EFlowTrack_1,  RecoedPrimaryTracks )")
            .Define("n_SecondaryTracks", "ReconstructedParticle2Track::getTK_n( SecondaryTracks )" )
            
            # which of the tracks are primary according to the reco algprithm
            .Define("IsPrimary_based_on_reco",  "VertexFitterSimple::IsPrimary_forTracks( EFlowTrack_1,  RecoedPrimaryTracks )")
            
            # jet clustering (ee-kt) before reconstructing SVs in event
            # build psedo-jets with the Reconstructed final particles
            .Define("pseudo_jets",    "JetClusteringUtils::set_pseudoJets_xyzm(RP_px, RP_py, RP_pz, RP_m)")
            # run jet clustering with all reco particles. ee_kt_algorithm, exclusive clustering, exactly 2 jets, E-scheme
            .Define("FCCAnalysesJets_ee_kt", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pseudo_jets)")
            # get the jets out of the structure
            .Define("jets_ee_kt", "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_kt)")
            # get the jet constituents out of the structure
            .Define("jetconstituents", "JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_kt)")
            # find SVs in jets
            .Define("SV_jet", "VertexFinderLCFIPlus::get_SV_jets(ReconstructedParticles, EFlowTrack_1, PV, IsPrimary_based_on_reco, jets_ee_kt, jetconstituents)")
            # find V0s in jets
            .Define("V0", "VertexFinderLCFIPlus::get_V0s_jet(ReconstructedParticles, EFlowTrack_1, IsPrimary_based_on_reco, jets_ee_kt, jetconstituents, PV)")
            # separate by jetsh. none of this will be needed once get_V0_jet is of the form Rvec<RVec<V0>>
            .Define("V0_jet", "VertexingUtils::get_svInJets(V0.vtx, V0.nSV_jet)")


            ###############
            #Some hacky reshaping of the constis 
            .Define("RPj_hasTRK",     "JetClusteringUtils::reshape2jet(RP_hasTRK, jetconstituents)")
            #.Define("split_indices", "ReconstructedParticle::index_splitter(ReconstructedParticles, RP_hasTRK)")
            .Define("split_indices",  "ReconstructedParticle::index_splitter(RP_hasTRK)")
            .Define("split_jetconstituents", "ReconstructedParticle::index_converter(jetconstituents, split_indices)")
            .Define("charged_constis","ReconstructedParticle::sel_template(true)(RPj_hasTRK, split_jetconstituents)")
            .Define("neutral_constis","ReconstructedParticle::sel_template(false)(RPj_hasTRK, split_jetconstituents)")
            
            #Some hacky reshaping of the event level RP ids
            .Define("RPj_isMuon",     "JetClusteringUtils::reshape2jet(RP_isMuon, jetconstituents)")
            .Define("RPj_isElectron", "JetClusteringUtils::reshape2jet(RP_isElectron, jetconstituents)")
            .Define("RPj_isPhoton",   "JetClusteringUtils::reshape2jet(RP_isPhoton, jetconstituents)")
                
            
            #Get jet-level vars
            .Define("jets_nRP_charged","JetClusteringUtils::get_nConstituents(charged_constis)")
            .Define("jets_nRP_neutral","JetClusteringUtils::get_nConstituents(neutral_constis)")
            .Define("jets_p",          "JetClusteringUtils::get_p(jets_ee_kt)")
            .Define("jets_px",         "JetClusteringUtils::get_px(jets_ee_kt)")
            .Define("jets_py",         "JetClusteringUtils::get_py(jets_ee_kt)")
            .Define("jets_pz",         "JetClusteringUtils::get_pz(jets_ee_kt)")
            .Define("jets_theta",      "JetClusteringUtils::get_theta(jets_ee_kt)")
            .Define("jets_phi",        "JetClusteringUtils::get_phi(jets_ee_kt)")
            .Define("jets_m",          "JetClusteringUtils::get_m(jets_ee_kt)")
            .Define("jets_e",          "JetClusteringUtils::get_e(jets_ee_kt)")

            #.Define("n_jets",          "JetClusteringUtils::get_n_jets(jets_ee_kt)")
            .Define("n_jets",          "jets_ee_kt.size()")
            
            .Define("jets_pt",         "JetClusteringUtils::get_pt(jets_ee_kt)")
            .Define("jets_eta",        "JetClusteringUtils::get_eta(jets_ee_kt)")
            
            .Define("jets_ghostFlavour", "JetTaggingUtils::get_flavour(Particle, Particle1, FCCAnalysesJets_ee_kt, pseudo_jets)")
            .Define("jets_flavour", "ROOT::VecOps::RVec<float> abs_flav; for(auto& flav : jets_ghostFlavour){if(flav==21){abs_flav.push_back(0); continue;};abs_flav.push_back(std::abs(flav));}; return abs_flav")

            #.Define("jets_flavour",    "JetTaggingUtils::get_flavour(jets_ee_kt, Particle)")
            
            #.Define("jets_flavour_placeholder","JetTaggingUtils::get_flavour(jets_ee_kt, Particle)")
            .Define("jets_onehot", "ReconstructedParticle::one_hot_encode(jets_flavour)")
            .Define("jets_isBbar", "jets_onehot[0]")
            .Define("jets_isCbar", "jets_onehot[1]")
            .Define("jets_isSbar", "jets_onehot[2]")
            .Define("jets_isUbar", "jets_onehot[3]")
            .Define("jets_isDbar", "jets_onehot[4]")
            .Define("jets_isD",    "jets_onehot[6]")
            .Define("jets_isU",    "jets_onehot[7]")
            .Define("jets_isS",    "jets_onehot[8]")
            .Define("jets_isC",    "jets_onehot[9]")
            .Define("jets_isB",    "jets_onehot[10]")
            .Define("jets_isUndefined","jets_onehot[5]")
            
            #.Define("jets_isBbar_placeholder", "jets_onehot[0]")
            #.Define("jets_isCbar_placeholder", "jets_onehot[1]")
            #.Define("jets_isSbar_placeholder", "jets_onehot[2]")
            #.Define("jets_isUbar_placeholder", "jets_onehot[3]")
            #.Define("jets_isDbar_placeholder", "jets_onehot[4]")
            #.Define("jets_isD_placeholder",    "jets_onehot[6]")
            #.Define("jets_isU_placeholder",    "jets_onehot[7]")
            #.Define("jets_isS_placeholder",    "jets_onehot[8]")
            #.Define("jets_isC_placeholder",    "jets_onehot[9]")
            #.Define("jets_isB_placeholder",    "jets_onehot[10]")
            #.Define("jets_isUndefined_placeholder","jets_onehot[5]")
               
            #.Define("jets_isBbar",     "int jets_isBbar = jets_isBbar_placeholder[0]; return jets_isBbar;")
            #.Define("jets_isCbar",     "int jets_isCbar = jets_isCbar_placeholder[0]; return jets_isCbar;")
            #.Define("jets_isSbar",     "int jets_isSbar = jets_isSbar_placeholder[0]; return jets_isSbar;")
            #.Define("jets_isUbar",     "int jets_isUbar = jets_isUbar_placeholder[0]; return jets_isUbar;")
            #.Define("jets_isDbar",     "int jets_isDbar = jets_isDbar_placeholder[0]; return jets_isDbar;")
            #.Define("jets_isD",        "int jets_isD = jets_isD_placeholder[0]; return jets_isD;")
            #.Define("jets_isU",        "int jets_isU = jets_isU_placeholder[0]; return jets_isU;")
            #.Define("jets_isS",        "int jets_isS = jets_isS_placeholder[0]; return jets_isS;")
            #.Define("jets_isC",        "int jets_isC = jets_isC_placeholder[0]; return jets_isC;")
            #.Define("jets_isB",        "int jets_isB = jets_isB_placeholder[0]; return jets_isB;")
            #.Define("jets_isUndefined","int jets_isUndefined = jets_isUndefined_placeholder[0]; return jets_isUndefined;") #is this undefined or u+d?
            
            
            #.Define("jets_ee_kt_e",          "JetClusteringUtils::get_e(jets_ee_kt)")
            #.Define("jets_ee_kt_px",         "JetClusteringUtils::get_px(jets_ee_kt)")
            #.Define("jets_ee_kt_py",         "JetClusteringUtils::get_py(jets_ee_kt)")
            #.Define("jets_ee_kt_pz",         "JetClusteringUtils::get_pz(jets_ee_kt)")
            #.Define("jets_ee_kt_flavour",    "JetTaggingUtils::get_flavour(jets_ee_kt, Particle)")
            #.Define("jets_ee_kt_flavour_anti","JetTaggingUtils::get_flavour_anti(jets_ee_kt, Particle)")
            #.Define("jets_theta_placeholder", "JetClusteringUtils::get_theta(jets_ee_kt)")
            #.Define("jets_phi_placeholder",   "JetClusteringUtils::get_phi(jets_ee_kt)")
            #.Define("jets_p_placeholder",     "JetClusteringUtils::get_p(jets_ee_kt)")


               
            #Reshape event-level charged RP to jet-level charged RP
            .Define("RPj_charged_p",    "JetClusteringUtils::reshape2jet(RP_charged_p, charged_constis)")
            .Define("RPj_charged_theta","JetClusteringUtils::reshape2jet(RP_charged_theta, charged_constis)")
            .Define("RPj_charged_phi",  "JetClusteringUtils::reshape2jet(RP_charged_phi, charged_constis)")
            
            #.Define("RPj_charged_theta_placeholder","JetClusteringUtils::reshape2jet(RP_charged_theta, charged_constis)")
            #.Define("RPj_charged_phi_placeholder",  "JetClusteringUtils::reshape2jet(RP_charged_phi, charged_constis)")
            #.Define("RPj_charged_p_placeholder",    "JetClusteringUtils::reshape2jet(RP_charged_p, charged_constis)")
            
            .Define("RPj_charged_mass",      "JetClusteringUtils::reshape2jet(RP_charged_mass, charged_constis)")
            .Define("RPj_charged_Z0",        "JetClusteringUtils::reshape2jet(RP_charged_Z0, charged_constis)")
            .Define("RPj_charged_D0",        "JetClusteringUtils::reshape2jet(RP_charged_D0, charged_constis)")
            .Define("RPj_charged_Z0_sig",    "JetClusteringUtils::reshape2jet(RP_charged_Z0_sig, charged_constis)")
            .Define("RPj_charged_D0_sig",    "JetClusteringUtils::reshape2jet(RP_charged_D0_sig, charged_constis)")
            .Define("RPj_charged_dTheta",    "JetClusteringUtils::get_dTheta(jets_theta, RPj_charged_theta)")
            .Define("RPj_charged_dPhi",      "JetClusteringUtils::get_dPhi(jets_phi, RPj_charged_phi)")
            .Define("RPj_charged_pRel",      "JetClusteringUtils::get_pRel(jets_p, RPj_charged_p)")
            .Define("RPj_charged_isMuon",    "ReconstructedParticle::sel_template(true)(RPj_hasTRK, RPj_isMuon)")
            .Define("RPj_charged_isElectron","ReconstructedParticle::sel_template(true)(RPj_hasTRK, RPj_isElectron)")
            
               
            #Reshape event-level neutral RP to jet-level neutral RP
            .Define("RPj_neutral_p", "JetClusteringUtils::reshape2jet(RP_neutral_p, neutral_constis)")
               
            #.Define("RPj_neutral_p_placeholder", "JetClusteringUtils::reshape2jet(RP_neutral_p, neutral_constis)")
            
            .Define("RPj_neutral_pRel",     "JetClusteringUtils::get_pRel(jets_p, RPj_neutral_p)")
            .Define("RPj_neutral_isPhoton", "ReconstructedParticle::sel_template(false)(RPj_hasTRK, RPj_isPhoton)")
            
            #Aliasing to match Alexandre's file
            .Define("isB", "std::vector<int> isB(jets_isB.begin(), jets_isB.end()); return isB;")
            .Define("isC", "std::vector<int> isC(jets_isC.begin(), jets_isC.end()); return isC;")
            .Define("isUD", "std::vector<int> isUD; for(int i=0; i<jets_isU.size(); ++i){isUD.push_back(jets_isU[i]+jets_isD[i]);}; return isUD;")
            .Define("isS", "std::vector<int> isS(jets_isS.begin(), jets_isS.end()); return isS;")
            #.Define("isUndefined", "std::vector<int> isUndefined(jets_isB.size(), 0); for(int i=0; i<isUndefined.size(); ++i){isUndefined[i] = (jets_isB[i]==0);}; return isUndefined;")
            #.Define("isUndefined", "std::vector<int> isUndefined(jets_isC.size(), 0); for(int i=0; i<isUndefined.size(); ++i){isUndefined[i] = (jets_isC[i]==0);}; return isUndefined;")
            .Define("isUndefined", "std::vector<int> isUndefined(jets_isU.size(), 0); for(int i=0; i<isUndefined.size(); ++i){isUndefined[i] = (jets_isU[i]==0)&&(jets_isD[i]==0)&&(jets_isS[i]==0);}; return isUndefined;")
            #.Define("isB", "int isB = jets_isB_placeholder[0]; return isB;")
            #.Define("isC", "int isC = jets_isC_placeholder[0]; return isC;")
            #.Define("isUD", "int isUD = jets_isU_placeholder[0]+jets_isD_placeholder[0]; return isUD;")
            #.Define("isS", "int isS = jets_isS_placeholder[0]; return isS;")
            #.Define("isUndefined", "int isUndefined = 0; return isUndefined;")
            .Define("jet_pt",                   "JetClusteringUtils::get_pt(jets_ee_kt)")
            .Define("jet_eta",                  "JetClusteringUtils::get_eta(jets_ee_kt)")
            .Define("n_Cpfcand",                "std::vector<int> n_Cpfcand(jets_nRP_charged.begin(), jets_nRP_charged.end()); return n_Cpfcand;")
            #.Define("n_Cpfcand",                "int n_Cpfcand = jets_nRP_charged; return n_Cpfcand;")
            .Define("nCpfcand",                 "jets_nRP_charged;")
            .Define("n_Npfcand",                "std::vector<int> n_Npfcand(jets_nRP_neutral.begin(), jets_nRP_neutral.end()); return n_Npfcand;")
            #.Define("n_Npfcand",                "int n_Npfcand = jets_nRP_neutral); return n_Npfcand;")
            .Define("nNpfcand",                 "jets_nRP_neutral")

            # what's going on? from here -
            .Define("Cpfcan_BtagPf_trackEtaRel","RPj_charged_p")
            .Define("Cpfcan_BtagPf_trackPtRel", "RPj_charged_theta")
            .Define("Cpfcan_BtagPf_trackDeltaR","RPj_charged_phi")
            .Define("Cpfcan_quality",           "RPj_charged_mass")
            .Define("Npfcan_ptrel",             "RPj_neutral_p")
            # - to here
            .Define("Npfcan_isGamma",           "RPj_neutral_isPhoton")
            .Define("isU",                      "std::vector<int> isU(jets_isU.begin(), jets_isU.end()); return isU;")
            .Define("isD",                      "std::vector<int> isD(jets_isD.begin(), jets_isD.end()); return isD;")
            #.Define("isU",                      "int isU = jets_isU[0]; return isU;")
            #.Define("isD",                      "int isD = jets_isD[0]; return isD;")

            .Define("Z_flavour", "JetTaggingUtils::get_Z_flavour(jets_ee_kt, Particle)")

            ###############



            
            # JET VARIABLES
            .Define("n_sv",       "VertexingUtils::get_n_SV_jets(SV_jet)") # no of SVs per jet

            # SV VARIABLES
            #.Define("sv_mass1",   "myUtils::get_Vertex_mass(SV.vtx, ReconstructedParticles)") # SV mass (first way)
            .Define("sv_mass",    "VertexingUtils::get_invM(SV_jet)") # SV mass (second way)
            #.Define("sv_p4",      "VertexingUtils::get_p4_SV(SV_jet)") # SV momentum (4 vector)
            .Define("sv_p",       "VertexingUtils::get_pMag_SV(SV_jet)") # SV momentum (magnitude)
            .Define("sv_ntracks", "VertexingUtils::get_VertexNtrk(SV_jet)") # SV daughters (no of tracks)
            .Define("sv_chi2",    "VertexingUtils::get_chi2_SV(SV_jet)") # SV chi2 (not normalised)
            .Define("sv_normchi2","VertexingUtils::get_norm_chi2_SV(SV_jet)") # SV chi2 (normalised)
            .Define("sv_ndf",     "VertexingUtils::get_nDOF_SV(SV_jet)") # SV no of DOF
            .Define("sv_theta",   "VertexingUtils::get_theta_SV(SV_jet)") # SV polar angle (theta)
            .Define("sv_phi",     "VertexingUtils::get_phi_SV(SV_jet)") # SV azimuthal angle (phi)
            .Define("sv_thetarel","VertexingUtils::get_relTheta_SV(SV_jet, jets_ee_kt)") # SV polar angle wrt jet
            .Define("sv_phirel",  "VertexingUtils::get_relPhi_SV(SV_jet, jets_ee_kt)") # SV azimuthal angle wrt jet
            .Define("sv_costhetasvpv","VertexingUtils::get_pointingangle_SV(SV_jet, PV)") # SV pointing angle
            .Define("sv_dxy",     "VertexingUtils::get_dxy_SV(SV_jet, PV)") # SV distance from PV (in xy plane)
            .Define("sv_d3d",     "VertexingUtils::get_d3d_SV(SV_jet, PV)") # SV distance from PV (in 3D)

            # V0 VARIABLES
            .Define("v0_pid",     "VertexingUtils::get_pdg_V0(V0.pdgAbs, V0.nSV_jet)") # V0 pdg id
            .Define("v0_mass",    "VertexingUtils::get_invM_V0(V0.invM, V0.nSV_jet)") # V0 invariant mass
            #.Define("v0_p4",      "VertexingUtils::get_p4_SV(V0_jet)") # V0 momentum (4 vector)
            .Define("v0_p",       "VertexingUtils::get_pMag_SV(V0_jet)") # V0 momentum (magnitude)
            .Define("v0_ntracks", "VertexingUtils::get_VertexNtrk(V0_jet)") # V0 daughters (no of tracks)
            .Define("v0_chi2",    "VertexingUtils::get_chi2_SV(V0_jet)") # V0 chi2 (not normalised)
            .Define("v0_normchi2","VertexingUtils::get_norm_chi2_SV(V0_jet)") # V0 chi2 (normalised but same as above)
            .Define("v0_ndf",     "VertexingUtils::get_nDOF_SV(V0_jet)") # V0 no of DOF (always 1)
            .Define("v0_theta",   "VertexingUtils::get_theta_SV(V0_jet)") # V0 polar angle (theta)
            .Define("v0_phi",     "VertexingUtils::get_phi_SV(V0_jet)") # V0 azimuthal angle (phi)
            .Define("v0_thetarel","VertexingUtils::get_relTheta_SV(V0_jet, jets_ee_kt)") # V0 polar angle wrt jets
            .Define("v0_phirel",  "VertexingUtils::get_relPhi_SV(V0_jet, jets_ee_kt)") # V0 azimuthal angle wrt jets
            .Define("v0_costhetasvpv","VertexingUtils::get_pointingangle_SV(V0_jet, PV)") # V0 pointing angle
            .Define("v0_dxy",     "VertexingUtils::get_dxy_SV(V0_jet, PV)") # V0 distance from PV (in xy plane)
            .Define("v0_d3d",     "VertexingUtils::get_d3d_SV(V0_jet, PV)") # V0 distance from PV (in 3D)

            
        )
        return df2

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [

            # charged PF
            "RPj_charged_p",
            "RPj_charged_theta",
            "RPj_charged_phi",
            "RPj_charged_mass",
            "RPj_charged_Z0",
            "RPj_charged_D0",
            "RPj_charged_Z0_sig",
            "RPj_charged_D0_sig",
            "RPj_charged_dTheta",
            "RPj_charged_dPhi",
            "RPj_charged_pRel",
            "RPj_charged_isMuon",
            "RPj_charged_isElectron",
            
            # neutral PF
            "RPj_neutral_p",
            "RPj_neutral_pRel",
            "RPj_neutral_isPhoton",            

            # jet variables
            "n_sv",

            "jets_nRP_charged",
            "jets_nRP_neutral",
            "jets_p",
            "jets_px",
            "jets_py",
            "jets_pz",
            "jets_theta",
            "jets_phi",
            "jets_m",
            "jets_e",

            "n_jets",
            
            "jets_pt",
            "jets_eta",
            
            "jets_isBbar",
            "jets_isCbar",
            "jets_isSbar",
            "jets_isUbar",
            "jets_isDbar",
            "jets_isD",
            "jets_isU",
            "jets_isS",
            "jets_isC",
            "jets_isB",
            "jets_isUndefined",
            
            "isB",
            "isC",
            "isUD",
            "isS",
            "isUndefined",
            "jet_pt",
            "jet_eta",
            "n_Cpfcand",
            "nCpfcand",
            "Cpfcan_BtagPf_trackEtaRel",
            "Cpfcan_BtagPf_trackPtRel",
            "Cpfcan_BtagPf_trackDeltaR",
            "Cpfcan_quality",
            "n_Npfcand",
            "nNpfcand",
            "Npfcan_ptrel",
            "Npfcan_isGamma",
            "isU",
            "isD",

            "Z_flavour",

            # sv variables
            #"sv_mass1",
            "sv_mass",
            #"sv_p4",
            "sv_p",
            "sv_ntracks",
            "sv_chi2",
            "sv_normchi2",
            "sv_ndf",
            "sv_theta",
            "sv_phi",
            "sv_thetarel",
            "sv_phirel",
            "sv_costhetasvpv",
            "sv_dxy",
            "sv_d3d",

            # v0 variables
            "v0_pid",
            "v0_mass",
            #"v0_p4",
            "v0_p",
            "v0_ntracks",
            "v0_chi2",
            "v0_normchi2",
            "v0_ndf",
            "v0_theta",
            "v0_phi",
            "v0_thetarel",
            "v0_phirel",
            "v0_costhetasvpv",
            "v0_dxy",
            "v0_d3d",

            ## ENDS HERE CURRENTLY

        ]
        return branchList
