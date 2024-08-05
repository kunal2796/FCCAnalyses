# testFile="root://eospublic.cern.ch//eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/p8_ee_Zbb_ecm91/events_066726720.root"

#Mandatory: List of processes  
processList = {
    # 'p8_ee_Zuds_ecm91':{'fraction':0.0001}, #Run 0.01% statistics in one output file named <outputDir>/p8_ee_Zuds_ecm91.root
    # 'p8_ee_Zcc_ecm91':{'fraction':0.0001}, #Run 0.01% statistics in one output file named <outputDir>/p8_ee_Zcc_ecm91.root
    # 'p8_ee_Zbb_ecm91':{'fraction':0.0001, 'output':'p8_ee_Zbb_ecm91_SV_100K'},
    # 'p8_ee_Zbb_ecm91':{'fraction':0.00001, 'output':'test'},
    'p8_ee_Zbb_ecm91':{'fraction':0.0001, 'output':'p8_ee_Zbb_ecm91'},
    'p8_ee_Zcc_ecm91':{'fraction':0.0001, 'output':'p8_ee_Zcc_ecm91'},
    # 'p8_ee_Zss_ecm91':{'fraction':0.0001, 'output':'p8_ee_Zss_ecm91'},
    # 'p8_ee_Zud_ecm91':{'fraction':0.0004, 'output':'p8_ee_Zud_ecm91'},
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/winter2023/IDEA/"

#Optional: output directory, default is local running directory
outputDir   = "/afs/cern.ch/work/k/kgautam/private/FCCAnalyses/outputs/FCCee/Zqq/asymmetries"

#Optional: ncpus, default is 4
nCPUS       = 8

#Optional running on HTCondor, default is False
# runBatch    = True

#Optional batch queue name when running on HTCondor, default is workday  
#batchQueue = "microcentury"
# batchQueue = "longlunch"

#Optional computing account when running on HTCondor, default is group_u_FCC.local_gen
# compGroup = "group_u_FCC.local_gen"

#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df2 = (df

               # alias
               .Alias("Particle0", "Particle#0.index") 
               .Alias("Particle1", "Particle#1.index") 
               .Alias("Jet3","Jet#3.index")
               .Alias("Muon0","Muon#0.index")
               .Alias("Electron0","Electron#0.index")
               .Alias("Photon0","Photon#0.index")
               .Alias("MCRecoAssociations0","MCRecoAssociations#0.index")
               .Alias("MCRecoAssociations1","MCRecoAssociations#1.index")
               .Define("Particle1_indices", "Particle1")

               # Reco'd particles
               .Define("RP_px",        "ReconstructedParticle::get_px(ReconstructedParticles)")
               .Define("RP_py",        "ReconstructedParticle::get_py(ReconstructedParticles)")
               .Define("RP_pz",        "ReconstructedParticle::get_pz(ReconstructedParticles)")               
               .Define("RP_e",         "ReconstructedParticle::get_e(ReconstructedParticles)")
               .Define("RP_m",         "ReconstructedParticle::get_mass(ReconstructedParticles)")
               .Define("RP_charge",    "ReconstructedParticle::get_charge(ReconstructedParticles)")
               .Define('RP_pdg',       "ReconstructedParticle2MC::getRP2MC_pdg(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")
               #.Define("RP_hasTRK",    "ReconstructedParticle2Track::hasTRK(ReconstructedParticles)")
               .Define("RP_hasTRK",    "ROOT::VecOps::RVec<int> result; for (auto& charge : RP_charge){if(std::abs(charge)>0) result.push_back(1); else result.push_back(0);}return result;")
               .Define("RP_charged",   "ReconstructedParticle::sel_tag(1)(RP_hasTRK, ReconstructedParticles)")
               .Define("RP_neutral",   "ReconstructedParticle::sel_tag(0)(RP_hasTRK, ReconstructedParticles)")

               # # MC PID
               # .Define("RP_PID",       "ReconstructedParticle::get_PID(MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle)")
               # #            
               # .Define("RP_isMuon",    "ReconstructedParticle::is_PID(13, RP_PID)")
               # .Define("RP_isElectron","ReconstructedParticle::is_PID(11, RP_PID)")
               # .Define("RP_isPion",    "ReconstructedParticle::is_PID(211, RP_PID)")
               # .Define("RP_isKaon",    "ReconstructedParticle::is_PID(321, RP_PID)")
               # .Define("RP_isProton",  "ReconstructedParticle::is_PID(2212, RP_PID)")
               # .Define("RP_isPhoton",  "ReconstructedParticle::is_PID(22, RP_PID)")
               # .Define("RP_isKlong",   "ReconstructedParticle::is_PID(130, RP_PID)")
               # .Define("RP_isNeutron", "ReconstructedParticle::is_PID(2112, RP_PID)")
               # .Define("RP_isNPion",   "ReconstructedParticle::is_PID(111, RP_PID)")

               # charged PF vars 
               .Define("RP_charged_charge",    "ReconstructedParticle::get_charge(RP_charged)")
               .Define("RP_charged_p",         "ReconstructedParticle::get_p(RP_charged)")
               #.Define("RP_charged_p_template","ReconstructedParticle::get_p(RP_template)")
               .Define("RP_charged_theta",     "ReconstructedParticle::get_theta(RP_charged)")
               .Define("RP_charged_phi",       "ReconstructedParticle::get_phi(RP_charged)")
               .Define("RP_charged_e",         "ReconstructedParticle::get_e(RP_charged)")
               .Define("RP_charged_mass",      "ReconstructedParticle::get_mass(RP_charged)") #CAUTION: true mass of RP or assuming pion?
               .Define("RP_charged_Z0",        "ReconstructedParticle2Track::getRP2TRK_Z0(RP_charged, EFlowTrack_1)")
               .Define("RP_charged_D0",        "ReconstructedParticle2Track::getRP2TRK_D0(RP_charged, EFlowTrack_1)")
               .Define("RP_charged_Z0_sig",    "ReconstructedParticle2Track::getRP2TRK_Z0_sig(RP_charged, EFlowTrack_1)")
               .Define("RP_charged_D0_sig",    "ReconstructedParticle2Track::getRP2TRK_D0_sig(RP_charged, EFlowTrack_1)")
               .Define("RP_charged_Curv",      "ReconstructedParticle2Track::getRP2TRK_omega(RP_charged, EFlowTrack_1)")

               # neutral PF vars 
               .Define("RP_neutral_p",    "ReconstructedParticle::get_p(RP_neutral)")
               .Define("RP_neutral_theta","ReconstructedParticle::get_theta(RP_neutral)")
               .Define("RP_neutral_phi",  "ReconstructedParticle::get_phi(RP_neutral)")
               .Define("RP_neutral_e",    "ReconstructedParticle::get_e(RP_neutral)")
               .Define("RP_neutral_mass", "ReconstructedParticle::get_mass(RP_neutral)") #CAUTION: even more strange for neutrals

               # PID vars
               .Define("RP_dNdx","ReconstructedTrack::tracks_dNdx(EFlowTrack_1, EFlowTrack_1, EFlowTrack, EFlowTrack_2)")
               .Define("RP_tof", "ReconstructedTrack::tracks_TOF(EFlowTrack_1, EFlowTrack_1, EFlowTrack, TrackerHits)")
               
               # number of tracks
               .Define("n_tracks","ReconstructedParticle2Track::getTK_n(EFlowTrack_1)")

               # PV and SV reconstruction and Jet Clustering
               # determime the primary (and secondary) tracks without using the MC-matching:
	       # Select the tracks that are reconstructed as primaries
               .Define("RecoedPrimaryTracks",  "VertexFitterSimple::get_PrimaryTracks( EFlowTrack_1, true, 4.5, 20e-3, 300, 0., 0., 0.)")
               .Define("n_RecoedPrimaryTracks","ReconstructedParticle2Track::getTK_n( RecoedPrimaryTracks )")
               # the final primary vertex :
               .Define("PV",  "VertexFitterSimple::VertexFitter_Tk ( 1, RecoedPrimaryTracks, true, 4.5, 20e-3, 300) ")
               .Define("PrimaryVertex",        "VertexingUtils::get_VertexData( PV )")
               .Define("PrimaryVertexPosition","PV.vertex.position")
               # the secondary tracks
               .Define("SecondaryTracks",    "VertexFitterSimple::get_NonPrimaryTracks( EFlowTrack_1,  RecoedPrimaryTracks )")
               .Define("n_SecondaryTracks",  "ReconstructedParticle2Track::getTK_n( SecondaryTracks )" )
               # which of the tracks are primary according to the reco algprithm (boolean)
               .Define("IsPrimary_based_on_reco",  "VertexFitterSimple::IsPrimary_forTracks( EFlowTrack_1,  RecoedPrimaryTracks )")

               # jet clustering (ee-kt) before reconstructing SVs in event
               #build psedo-jets with the Reconstructed final particles
               .Define("pseudo_jets", "JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")
               #run jet clustering with all reco particles. ee_kt_algorithm, exclusive clustering, exactly 2 jets, E-scheme
               .Define("FCCAnalysesJets_ee_kt", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pseudo_jets)")
               #get the jets out of the structure
               .Define("jets_ee_kt", "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_kt)")
               #get the jet constituents out of the structure
	       .Define("jetconstituents", "JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_kt)")

	       # finding SVs in jets 
               .Define("SV_jet", "VertexFinderLCFIPlus::get_SV_jets(ReconstructedParticles, EFlowTrack_1, PV, IsPrimary_based_on_reco, jets_ee_kt, jetconstituents)")
               #find V0s in jets
               .Define("V0", "VertexFinderLCFIPlus::get_V0s_jet(ReconstructedParticles, EFlowTrack_1, IsPrimary_based_on_reco, jets_ee_kt, jetconstituents, PV)")
               #separate by jetsh. none of this will be needed once get_V0_jet is of the form Rvec<RVec<V0>>
               .Define("V0_jet", "VertexingUtils::get_svInJets(V0.vtx, V0.nSV_jet)")
               #separate tracks by jets
               .Define("Tracks", "VertexingUtils::get_tracksInJets(ReconstructedParticles, EFlowTrack_1, jets_ee_kt, jetconstituents)")

               # jet-level variables
               # .Define("jets_nRP_charged","JetClusteringUtils::get_nConstituents(charged_constis)") #won't work due to other fns unavailable to get charged_constis
               # .Define("jets_nRP_neutral","JetClusteringUtils::get_nConstituents(neutral_constis)")
               .Define("jets_p",          "JetClusteringUtils::get_p(jets_ee_kt)")
               .Define("jets_px",         "JetClusteringUtils::get_px(jets_ee_kt)")
               .Define("jets_py",         "JetClusteringUtils::get_py(jets_ee_kt)")
               .Define("jets_pz",         "JetClusteringUtils::get_pz(jets_ee_kt)")
               .Define("jets_theta",      "JetClusteringUtils::get_theta(jets_ee_kt)")
               .Define("jets_phi",        "JetClusteringUtils::get_phi(jets_ee_kt)")
               .Define("jets_m",          "JetClusteringUtils::get_m(jets_ee_kt)")
               .Define("jets_e",          "JetClusteringUtils::get_e(jets_ee_kt)")
               .Define("jets_pt",         "JetClusteringUtils::get_pt(jets_ee_kt)")
               .Define("jets_eta",        "JetClusteringUtils::get_eta(jets_ee_kt)")
               .Define("n_jets",          "jets_ee_kt.size()")

               # jet flavour
               .Define("jets_flavour",  "JetTaggingUtils::get_flavour(jets_ee_kt, Particle)")
               .Define("Z_flavour",     "JetTaggingUtils::get_Zqq_flavour(jets_ee_kt, Particle)") # default signed flavour to be used (q-qbar separated)
               # .Redefine("jets_flavour","ROOT::VecOps::RVec<float> flav; for(int i=0; i<jets_flavour.size(); ++i){if(jets_flavour.at(i)==Z_flavour.at(i)) flav.push_back(jets_flavour.at(i)); else flav.push_back(0);}return flav;") # includes undefined flavour when q(qbar) is outside the cone

               # SV variables
               .Define("n_sv",       "VertexingUtils::get_n_SV_jets(SV_jet)") # no of SVs per jet
               .Define("sv_mass",    "VertexingUtils::get_invM(SV_jet)") # SV mass (second way)
               # .Define("sv_p4",      "VertexingUtils::get_p4_SV(SV_jet)") # SV momentum (4 vector)
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
               .Define("sv_trackIndx", "VertexingUtils::get_VertexRecoInd_SVjet(SV_jet, ReconstructedParticles)") #IMPORTANT AND TO BE FIXED

               # V0 variables
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
               
               ############### all reshaping redefining here #############
               # # Some hacky reshaping of the constis 
               # .Define("RPj_hasTRK",     "JetClusteringUtils::reshape2jet(RP_hasTRK, jetconstituents)")
               # .Define("split_indices",  "ReconstructedParticle::index_splitter(RP_hasTRK)")
               # .Define("split_jetconstituents", "ReconstructedParticle::index_converter(jetconstituents, split_indices)")
               # .Define("charged_constis","ReconstructedParticle::sel_template(true)(RPj_hasTRK, split_jetconstituents)")
               # .Define("neutral_constis","ReconstructedParticle::sel_template(false)(RPj_hasTRK, split_jetconstituents)")

               # #Some hacky reshaping of the event level RP ids
               # .Define("RPj_isMuon",     "JetClusteringUtils::reshape2jet(RP_isMuon, jetconstituents)")
               # .Define("RPj_isElectron", "JetClusteringUtils::reshape2jet(RP_isElectron, jetconstituents)")
               # .Define("RPj_isPion", "JetClusteringUtils::reshape2jet(RP_isPion, jetconstituents)")
               # .Define("RPj_is_Kaon", "JetClusteringUtils::reshape2jet(RP_isKaon, jetconstituents)")
               # .Define("RPj_isProton", "JetClusteringUtils::reshape2jet(RP_isProton, jetconstituents)")
               # .Define("RPj_isPhoton", "JetClusteringUtils::reshape2jet(RP_isPhoton, jetconstituents)")
               # .Define("RPj_isKlong", "JetClusteringUtils::reshape2jet(RP_isKlong, jetconstituents)")
               # .Define("RPj_isNeutron", "JetClusteringUtils::reshape2jet(RP_isNeutron, jetconstituents)")
               # .Define("RPj_isNPion", "JetClusteringUtils::reshape2jet(RP_isNPion, jetconstituents)")
               # .Define("RPj_PID", "JetClusteringUtils::reshape2jet(RP_PID, jetconstituents)")

               # # jet-level variables
               # .Define("jets_nRP_charged","JetClusteringUtils::get_nConstituents(charged_constis)")
               # .Define("jets_nRP_neutral","JetClusteringUtils::get_nConstituents(neutral_constis)")

               # # PID variables
               # .Define("RPj_dNdx","JetClusteringUtils::reshape2jet(RP_dNdx, jetconstituents)")
               # .Define("RPj_tof", "JetClusteringUtils::reshape2jet(RP_tof, jetconstituents)")

               # # kinematic variables of jet constituents
               # .Define("RPj_e",     "JetClusteringUtils::reshape2jet(RP_e, jetconstituents)")
               # .Define("RPj_px",    "JetClusteringUtils::reshape2jet(RP_px, jetconstituents)")
               # .Define("RPj_py",    "JetClusteringUtils::reshape2jet(RP_py, jetconstituents)")
               # .Define("RPj_pz",    "JetClusteringUtils::reshape2jet(RP_pz, jetconstituents)")
               # .Define("RPj_dAngle","JetClusteringUtils::get_dAngle(jets_px, jets_py, jets_pz, RPj_px, RPj_py, RPj_pz)")

               # # charged PF
               # .Define("RPj_charged_p",     "JetClusteringUtils::reshape2jet(RP_charged_p, charged_constis)")
               # .Define("RPj_charged_theta", "JetClusteringUtils::reshape2jet(RP_charged_theta, charged_constis)")
               # .Define("RPj_charged_phi",   "JetClusteringUtils::reshape2jet(RP_charged_phi, charged_constis)")
               # .Define("RPj_charged_e",     "JetClusteringUtils::reshape2jet(RP_charged_e, charged_constis)")
               # .Define("RPj_charged_mass",  "JetClusteringUtils::reshape2jet(RP_charged_mass, charged_constis)")
               # .Define("RPj_charged_charge","JetClusteringUtils::reshape2jet(RP_charged_charge, charged_constis)")
               # .Define("RPj_charged_Z0",    "JetClusteringUtils::reshape2jet(RP_charged_Z0, charged_constis)")
               # .Define("RPj_charged_D0",    "JetClusteringUtils::reshape2jet(RP_charged_D0, charged_constis)")
               # .Define("RPj_charged_Z0_sig","JetClusteringUtils::reshape2jet(RP_charged_Z0_sig, charged_constis)")
               # .Define("RPj_charged_D0_sig","JetClusteringUtils::reshape2jet(RP_charged_D0_sig, charged_constis)")
               # .Define("RPj_charged_Curv",  "JetClusteringUtils::reshape2jet(RP_charged_Curv, charged_constis)")
               # .Define("RPj_charged_pRel",  "JetClusteringUtils::get_pRel(jets_p, RPj_charged_p)")
               # .Define("RPj_charged_eRel",  "JetClusteringUtils::get_eRel(jets_p, RPj_charged_p)")
               # .Define("RPj_charged_dTheta","JetClusteringUtils::get_dTheta(jets_theta, RPj_charged_theta)")
               # .Define("RPj_charged_dPhi",  "JetClusteringUtils::get_dPhi(jets_phi, RPj_charged_phi)")
               # .Define("RPj_charged_dR",    "JetClusteringUtils::get_dR(RPj_charged_dTheta, RPj_charged_dPhi)")
               # .Define("RPj_charged_dAngle","ReconstructedParticle::sel_template(true)(RPj_hasTRK, RPj_dAngle)")
               # .Define("RPj_charged_isMuon",    "ReconstructedParticle::sel_template(true)(RPj_hasTRK, RPj_isMuon)")
               # .Define("RPj_charged_isElectron","ReconstructedParticle::sel_template(true)(RPj_hasTRK, RPj_isElectron)")
               # .Define("RPj_charged_isPion","ReconstructedParticle::sel_template(true)(RPj_hasTRK, RPj_isPion)")
               # .Define("RPj_charged_is_Kaon", "ReconstructedParticle::sel_template(true)(RPj_hasTRK, RPj_is_Kaon)")
               # .Define("RPj_charged_isProton","ReconstructedParticle::sel_template(true)(RPj_hasTRK, RPj_isProton)")
               # .Define("RPj_charged_PID", "ReconstructedParticle::sel_template(true)(RPj_hasTRK, RPj_PID)")
               
               # # neutral PF
               # .Define("RPj_neutral_p",     "JetClusteringUtils::reshape2jet(RP_neutral_p, neutral_constis)")
               # .Define("RPj_neutral_theta", "JetClusteringUtils::reshape2jet(RP_neutral_theta, neutral_constis)")
               # .Define("RPj_neutral_phi",   "JetClusteringUtils::reshape2jet(RP_neutral_phi, neutral_constis)")
               # .Define("RPj_neutral_e",     "JetClusteringUtils::reshape2jet(RP_neutral_e, neutral_constis)")
               # .Define("RPj_neutral_mass",  "JetClusteringUtils::reshape2jet(RP_neutral_mass, neutral_constis)")
               # .Define("RPj_neutral_pRel",  "JetClusteringUtils::get_pRel(jets_p, RPj_neutral_p)")
               # .Define("RPj_neutral_eRel",  "JetClusteringUtils::get_eRel(jets_p, RPj_neutral_p)")
               # .Define("RPj_neutral_dTheta","JetClusteringUtils::get_dTheta(jets_theta, RPj_neutral_theta)")
               # .Define("RPj_neutral_dPhi",  "JetClusteringUtils::get_dPhi(jets_phi, RPj_neutral_phi)")
               # .Define("RPj_neutral_dR",    "JetClusteringUtils::get_dR(RPj_neutral_dTheta, RPj_neutral_dPhi)")
               # .Define("RPj_neutral_dAngle","ReconstructedParticle::sel_template(false)(RPj_hasTRK, RPj_dAngle)")
               # .Define("RPj_neutral_isPhoton", "ReconstructedParticle::sel_template(false)(RPj_hasTRK, RPj_isPhoton)")
               # .Define("RPj_neutral_isKlong", "ReconstructedParticle::sel_template(false)(RPj_hasTRK, RPj_isKlong)")
               # .Define("RPj_neutral_isNeutron", "ReconstructedParticle::sel_template(false)(RPj_hasTRK, RPj_isNeutron)")
               # .Define("RPj_neutral_isNPion", "ReconstructedParticle::sel_template(false)(RPj_hasTRK, RPj_isNPion)")

               # .Define("n_Cpfcand","std::vector<int> n_Cpfcand(jets_nRP_charged.begin(), jets_nRP_charged.end()); return n_Cpfcand;")
               # .Define("nCpfcand", "jets_nRP_charged;") # I suppose the redundancy is for sanity checks
               # .Define("n_Npfcand","std::vector<int> n_Npfcand(jets_nRP_neutral.begin(), jets_nRP_neutral.end()); return n_Npfcand;")
               # .Define("nNpfcand", "jets_nRP_neutral")

               # #Redefinitions for consistency...
               # .Redefine("Z_flavour", "return std::vector<int>(Z_flavour.begin(), Z_flavour.end());")
               # .Redefine("jets_p",    "return std::vector<float>(jets_p.begin(), jets_p.end());")
               # .Redefine("jets_px",   "return std::vector<float>(jets_px.begin(), jets_px.end());")
               # .Redefine("jets_py",   "return std::vector<float>(jets_py.begin(), jets_py.end());")
               # .Redefine("jets_pz",   "return std::vector<float>(jets_pz.begin(), jets_pz.end());")
               # .Redefine("jets_theta","return std::vector<float>(jets_theta.begin(), jets_theta.end());")
               # .Redefine("jets_phi",  "return std::vector<float>(jets_phi.begin(), jets_phi.end());")
               # .Redefine("jets_m",    "return std::vector<float>(jets_m.begin(), jets_m.end());")
               # .Redefine("jets_e",    "return std::vector<float>(jets_e.begin(), jets_e.end());")
               # .Redefine("jets_nRP_charged", "return std::vector<float>(jets_nRP_charged.begin(), jets_nRP_charged.end());")
               # .Redefine("jets_nRP_neutral", "return std::vector<float>(jets_nRP_neutral.begin(), jets_nRP_neutral.end());")
               # #
               # .Redefine("sv_mass", "std::vector<std::vector<float>> result; for(auto& sv_mass_single : sv_mass){std::vector<float> tmp_res(sv_mass_single.begin(), sv_mass_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("sv_p", "std::vector<std::vector<float>> result; for(auto& sv_p_single : sv_p){std::vector<float> tmp_res(sv_p_single.begin(), sv_p_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("sv_ntracks", "std::vector<std::vector<int>> result; for(auto& sv_ntracks_single : sv_ntracks){std::vector<int> tmp_res(sv_ntracks_single.begin(), sv_ntracks_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("sv_chi2", "std::vector<std::vector<float>> result; for(auto& sv_chi2_single : sv_chi2){std::vector<float> tmp_res(sv_chi2_single.begin(), sv_chi2_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("sv_normchi2", "std::vector<std::vector<float>> result; for(auto& sv_normchi2_single : sv_normchi2){std::vector<float> tmp_res(sv_normchi2_single.begin(), sv_normchi2_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("sv_ndf", "std::vector<std::vector<int>> result; for(auto& sv_ndf_single : sv_ndf){std::vector<int> tmp_res(sv_ndf_single.begin(), sv_ndf_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("sv_theta", "std::vector<std::vector<float>> result; for(auto& sv_theta_single : sv_theta){std::vector<float> tmp_res(sv_theta_single.begin(), sv_theta_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("sv_phi", "std::vector<std::vector<float>> result; for(auto& sv_phi_single : sv_phi){std::vector<float> tmp_res(sv_phi_single.begin(), sv_phi_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("sv_thetarel", "std::vector<std::vector<float>> result; for(auto& sv_thetarel_single : sv_thetarel){std::vector<float> tmp_res(sv_thetarel_single.begin(), sv_thetarel_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("sv_phirel", "std::vector<std::vector<float>> result; for(auto& sv_phirel_single : sv_phirel){std::vector<float> tmp_res(sv_phirel_single.begin(), sv_phirel_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("sv_costhetasvpv", "std::vector<std::vector<float>> result; for(auto& sv_costhetasvpv_single : sv_costhetasvpv){std::vector<float> tmp_res(sv_costhetasvpv_single.begin(), sv_costhetasvpv_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("sv_dxy", "std::vector<std::vector<float>> result; for(auto& sv_dxy_single : sv_dxy){std::vector<float> tmp_res(sv_dxy_single.begin(), sv_dxy_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("sv_d3d", "std::vector<std::vector<float>> result; for(auto& sv_d3d_single : sv_d3d){std::vector<float> tmp_res(sv_d3d_single.begin(), sv_d3d_single.end()); result.push_back(tmp_res);} return result;")
               # #
               # .Redefine("v0_pid", "std::vector<std::vector<float>> result; for(auto& v0_pid_single : v0_pid){std::vector<float> tmp_res(v0_pid_single.begin(), v0_pid_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("v0_mass", "std::vector<std::vector<float>> result; for(auto& v0_mass_single : v0_mass){std::vector<float> tmp_res(v0_mass_single.begin(), v0_mass_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("v0_p", "std::vector<std::vector<float>> result; for(auto& v0_p_single : v0_p){std::vector<float> tmp_res(v0_p_single.begin(), v0_p_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("v0_ntracks", "std::vector<std::vector<int>> result; for(auto& v0_ntracks_single : v0_ntracks){std::vector<int> tmp_res(v0_ntracks_single.begin(), v0_ntracks_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("v0_chi2", "std::vector<std::vector<float>> result; for(auto& v0_chi2_single : v0_chi2){std::vector<float> tmp_res(v0_chi2_single.begin(), v0_chi2_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("v0_normchi2", "std::vector<std::vector<float>> result; for(auto& v0_normchi2_single : v0_normchi2){std::vector<float> tmp_res(v0_normchi2_single.begin(), v0_normchi2_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("v0_ndf", "std::vector<std::vector<int>> result; for(auto& v0_ndf_single : v0_ndf){std::vector<int> tmp_res(v0_ndf_single.begin(), v0_ndf_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("v0_theta", "std::vector<std::vector<float>> result; for(auto& v0_theta_single : v0_theta){std::vector<float> tmp_res(v0_theta_single.begin(), v0_theta_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("v0_phi", "std::vector<std::vector<float>> result; for(auto& v0_phi_single : v0_phi){std::vector<float> tmp_res(v0_phi_single.begin(), v0_phi_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("v0_thetarel", "std::vector<std::vector<float>> result; for(auto& v0_thetarel_single : v0_thetarel){std::vector<float> tmp_res(v0_thetarel_single.begin(), v0_thetarel_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("v0_phirel", "std::vector<std::vector<float>> result; for(auto& v0_phirel_single : v0_phirel){std::vector<float> tmp_res(v0_phirel_single.begin(), v0_phirel_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("v0_costhetasvpv", "std::vector<std::vector<float>> result; for(auto& v0_costhetasvpv_single : v0_costhetasvpv){std::vector<float> tmp_res(v0_costhetasvpv_single.begin(), v0_costhetasvpv_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("v0_dxy", "std::vector<std::vector<float>> result; for(auto& v0_dxy_single : v0_dxy){std::vector<float> tmp_res(v0_dxy_single.begin(), v0_dxy_single.end()); result.push_back(tmp_res);} return result;")
               # .Redefine("v0_d3d", "std::vector<std::vector<float>> result; for(auto& v0_d3d_single : v0_d3d){std::vector<float> tmp_res(v0_d3d_single.begin(), v0_d3d_single.end()); result.push_back(tmp_res);} return result;")

               
               
               )
        return df2

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [

            # jet flavour
            "Z_flavour",
            "jets_flavour",

            # jet-level variables
            "n_jets",
            "jets_p",
            "jets_px",
            "jets_py",
            "jets_pz",
            "jets_theta",
            "jets_phi",
            "jets_m",
            "jets_e",
            "jets_pt",
            "jets_eta",
            # "jets_nRP_charged",
            # "jets_nRP_neutral",
            "jetconstituents",

            # RP variables
            "RP_px",
            "RP_py",
            "RP_pz",
            "RP_e",
            "RP_m",
            "RP_charge",
            "RP_hasTRK",
            "RP_charged",
            "RP_neutral",

            # PID variables
            #"RPj_PID", #?
            # "RPj_dNdx",
            # "RPj_tof",
            "RP_pdg",
            "RP_dNdx",
            "RP_tof",

            # charged PF vars
            "n_tracks",
            "RP_charged_charge",
            "RP_charged_p",
            "RP_charged_theta",
            "RP_charged_phi",
            "RP_charged_e",
            "RP_charged_mass",
            "RP_charged_Z0",
            "RP_charged_D0",
            "RP_charged_Z0_sig",
            "RP_charged_D0_sig",
            "RP_charged_Curv",
            
            # neutral PF vars 
            "RP_neutral_p",
            "RP_neutral_theta",
            "RP_neutral_phi",
            "RP_neutral_e",
            "RP_neutral_mass",
            
            # # charged PF
            # "RPj_charged_p",
            # "RPj_charged_theta",
            # "RPj_charged_phi",
            # "RPj_charged_mass",
            # "RPj_charged_charge",
            # "RPj_charged_Z0",
            # "RPj_charged_D0",
            # "RPj_charged_Z0_sig",
            # "RPj_charged_D0_sig",
            # "RPj_charged_Curv",
            # "RPj_charged_pRel",
            # "RPj_charged_eRel",
            # "RPj_charged_dTheta",
            # "RPj_charged_dPhi",
            # "RPj_charged_dR",
            # "RPj_charged_dAngle",
            # "RPj_charged_isMuon",
            # "RPj_charged_isElectron",
            # "RPj_charged_isPion",
            # "RPj_charged_isKaon",
            # "RPj_charged_isProton",
            # "RPj_charged_PID",

            # # neutral PF
            # "RPj_neutral_p",
            # "RPj_neutral_theta",
            # "RPj_neutral_phi",
            # "RPj_neutral_mass",
            # "RPj_neutral_pRel",
            # "RPj_neutral_eRel",
            # "RPj_neutral_dTheta",
            # "RPj_neutral_dPhi",
            # "RPj_neutral_dR",
            # "RPj_neutral_dAngle",
            # "RPj_neutral_isPhoton",            
            # "RPj_neutral_isKlong",            
            # "RPj_neutral_isNeutron",            
            # "RPj_neutral_isNPion",

            # PV variables
            "PrimaryVertexPosition",
            "IsPrimary_based_on_reco",

            # SV variables
            #"nSV", # should add these. would require vec instead of vec of vec like for others. add while adding PID
            "n_sv",
            "sv_mass",
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
            "sv_trackIndx", #IMPORTANT AND TO BE FIXED
            
            #V0 variables
            #"nV0",
            "v0_pid",
            "v0_mass",
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
            
                ]
        return branchList
