#Mandatory: List of processes
processList = {
    'p8_ee_Zbb_ecm91':{'fraction':0.0001}, #Run 0.01% statistics in one output file named <outputDir>/p8_ee_Zbb_ecm91.root
    #'p8_ee_Zbb_ecm91':{'fraction':0.0001, 'output':'p8_ee_Zbb_ecm91_SV_100K'},
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/spring2021/IDEA/"

#Optional: output directory, default is local running directory
outputDir   = "outputs/FCCee/KG/"

#Optional: ncpus, default is 4
nCPUS       = 8

#Optional running on HTCondor, default is False
#runBatch    = False

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
            # global variables
            #.Define("")

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
            .Define("n_SecondaryTracks",  "ReconstructedParticle2Track::getTK_n( SecondaryTracks )" )
            
            # which of the tracks are primary according to the reco algprithm
            .Define("IsPrimary_based_on_reco",  "VertexFitterSimple::IsPrimary_forTracks( EFlowTrack_1,  RecoedPrimaryTracks )")
            
            # jet clustering (ee-kt) before reconstructing SVs in event
            .Define("RP_px",  "ReconstructedParticle::get_px(ReconstructedParticles)")
            .Define("RP_py",  "ReconstructedParticle::get_py(ReconstructedParticles)")
            .Define("RP_pz",  "ReconstructedParticle::get_pz(ReconstructedParticles)")
            .Define("RP_e",   "ReconstructedParticle::get_e(ReconstructedParticles)")
            # build psedo-jets with the Reconstructed final particles
            .Define("pseudo_jets", "JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")
            # run jet clustering with all reco particles. ee_kt_algorithm, exclusive clustering, exactly 2 jets, E-scheme
            .Define("FCCAnalysesJets_ee_kt", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pseudo_jets)")
            # get the jets out of the structure
            .Define("jets_ee_kt", "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_kt)")
            # get the jet constituents out of the structure
            .Define("jetconstituents_ee_kt", "JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_kt)")
            # find SVs in jets
            .Define("SV", "VertexFitterSimple::get_SV_jets(ReconstructedParticles, EFlowTrack_1, PV, IsPrimary_based_on_reco, jets_ee_kt, jetconstituents_ee_kt)")
            # find V0s in jets (yet to write the fn)
            .Define("V0", "VertexFitterSimple::get_V0_jet(ReconstructedParticles, EFlowTrack_1, IsPrimary_based_on_reco, jets_ee_kt, jetconstituents_ee_kt, PV)")

            # JET VARIABLES
            .Define("n_sv",       "SV.nSV_jet") # no of SVs per jet
            .Define("jet_p",      "JetClusteringUtils::get_p(jets_ee_kt)") # jet momentum
            .Define("jet_pt",     "JetClusteringUtils::get_pt(jets_ee_kt)") # jet transverse momentum
            .Define("jet_energy", "JetClusteringUtils::get_e(jets_ee_kt)") # jet energy
            .Define("jet_mass",   "JetClusteringUtils::get_m(jets_ee_kt)") # jet mass
            .Define("jet_theta",  "JetClusteringUtils::get_theta(jets_ee_kt)") # jet polar angle
            .Define("jet_eta",    "JetClusteringUtils::get_eta(jets_ee_kt)") # jet pseudo-rapidity
            .Define("jet_phi",    "JetClusteringUtils::get_phi(jets_ee_kt)") # jet azimuthal angle

            # SV VARIABLES
            .Define("sv_mass1",   "myUtils::get_Vertex_mass(SV.vtx, ReconstructedParticles)") # SV mass (first way)
            .Define("sv_mass",    "VertexingUtils::get_invM(SV.vtx)") # SV mass (second way)
            .Define("sv_p3",      "VertexingUtils::get_p_SV(SV.vtx)") # SV momentum (3 vector)
            .Define("sv_p",       "VertexingUtils::get_pMag_SV(SV.vtx)") # SV momentum (magnitude)
            .Define("sv_ntracks", "VertexingUtils::get_VertexNtrk(SV.vtx)") # SV daughters (no of tracks)
            .Define("sv_chi2",    "VertexingUtils::get_chi2_SV(SV)") # SV chi2 (not normalised)
            .Define("sv_normchi2","VertexingUtils::get_norm_chi2_SV(SV.vtx)") # SV chi2 (normalised)
            .Define("sv_ndf",     "VertexingUtils::get_nDOF_SV(SV.vtx)") # SV no of DOF
            .Define("sv_theta",   "VertexingUtils::get_theta_SV(SV.vtx)") # SV polar angle (theta)
            .Define("sv_phi",     "VertexingUtils::get_phi_SV(SV.vtx)") # SV azimuthal angle (phi)
            .Define("sv_thetarel","VertexingUtils::get_relTheta_SV(SV.vtx, SV.nSV_jet, jets_ee_kt)") # SV polar angle wrt jets
            .Define("sv_phirel",  "VertexingUtils::get_relPhi_SV(SV.vtx, SV.nSV_jet, jets_ee_kt)") # SV azimuthal angle wrt jets
            .Define("sv_costhetasvpv","VertexingUtils::get_pointingangle_SV(SV.vtx, PV)") # SV pointing angle
            .Define("sv_dxy",     "VertexingUtils::get_dxy_SV(SV.vtx, PV)") # SV distance from PV (in xy plane)
            .Define("sv_d3d",     "VertexingUtils::get_d3d_SV(SV.vtx, PV)") # SV distance from PV (in 3D)

            # V0 VARIABLES
            .Define("v0_pid",     "VertexingUtils::get_pdg_V0(V0)") # V0 pdg id
            .Define("v0_mass1",   "myUtils::get_Vertex_mass(V0.vtx, ReconstructedParticles)") # V0 mass (first way)
            .Define("v0_mass",    "VertexingUtils::get_invM_V0(V0)") # V0 mass (second way)
            .Define("v0_p3",      "VertexingUtils::get_p_SV(V0.vtx)") # V0 momentum (3 vector)
            .Define("v0_p",       "VertexingUtils::get_pMag_SV(V0.vtx)") # V0 momentum (magnitude)
            .Define("v0_ntracks", "VertexingUtils::get_VertexNtrk(V0.vtx)") # V0 daughters (no of tracks)
            .Define("v0_chi2",    "VertexingUtils::get_chi2_SV(V0)") # V0 chi2 (not normalised)
            .Define("v0_normchi2","VertexingUtils::get_norm_chi2_SV(V0.vtx)") # V0 chi2 (normalised but same as above)
            .Define("v0_ndf",     "VertexingUtils::get_nDOF_SV(V0.vtx)") # V0 no of DOF (always 1)
            .Define("v0_theta",   "VertexingUtils::get_theta_SV(V0.vtx)") # V0 polar angle (theta)
            .Define("v0_phi",     "VertexingUtils::get_phi_SV(V0.vtx)") # V0 azimuthal angle (phi)
            #.Define("v0_thetarel","VertexingUtils::get_relTheta_SV(V0.vtx, V0.nV0_jet, jets_ee_kt)") # V0 polar angle wrt jets
            #.Define("v0_phirel",  "VertexingUtils::get_relPhi_SV(V0.vtx, V0.nV0_jet, jets_ee_kt)") # V0 azimuthal angle wrt jets
            .Define("v0_costhetasvpv","VertexingUtils::get_pointingangle_SV(V0.vtx, PV)") # V0 pointing angle
            .Define("v0_dxy",     "VertexingUtils::get_dxy_SV(V0.vtx, PV)") # V0 distance from PV (in xy plane)
            .Define("v0_d3d",     "VertexingUtils::get_d3d_SV(V0.vtx, PV)") # V0 distance from PV (in 3D)

            # TRACK VARIABLES
            .Define("gtrack_dz",  "ReconstructedParticle2Track::getRP2TRK_Z0(ReconstructedParticles, EFlowTrack_1)") # longitudinal IP
            .Define("gtrack_dxy", "ReconstructedParticle2Track::getRP2TRK_D0(ReconstructedParticles, EFlowTrack_1)") # transverse IP
            .Define("gtrack_sdz", "ReconstructedParticle2Track::getRP2TRK_Z0_sig(ReconstructedParticles, EFlowTrack_1)") # longitudinal IP significance
            .Define("gtrack_sdxy","ReconstructedParticle2Track::getRP2TRK_D0_sig(ReconstructedParticles, EFlowTrack_1)") # transverse IP significance
            .Define("gtrack_phi", "ReconstructedParticle2Track::getRP2TRK_phi(ReconstructedParticles, EFlowTrack_1)") # azimuthal angle
            .Define("gtrack_theta","ReconstructedParticle2Track::getRP2TRK_theta(ReconstructedParticles, EFlowTrack_1)") # polar angle
            .Define("gtrack_p",   "ReconstructedParticle2Track::getRP2TRK_mom(ReconstructedParticles, EFlowTrack_1)") # momentum magnitude
            
        )
        return df2

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
            # jet variables
            "n_sv",
            "jet_p",
            "jet_pt",
            "jet_energy",
            "jet_mass",
            "jet_theta",
            "jet_eta",
            "jet_phi",

            # sv variables
            "sv_mass1",
            "sv_mass",
            "sv_p3",
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
            "v0_mass1",
            "v0_mass",
            "v0_p3",
            "v0_p",
            "v0_ntracks",
            "v0_chi2",
            "v0_normchi2",
            "v0_ndf",
            "v0_theta",
            "v0_phi",
            #"v0_thetarel",
            #"v0_phirel",
            "v0_costhetasvpv",
            "v0_dxy",
            "v0_d3d",

            # track variables
            "gtrack_dz",
            "gtrack_dxy",
            "gtrack_sdz",
            "gtrack_sdxy",
            "gtrack_phi",
            "gtrack_theta",
            "gtrack_p"
        ]
        return branchList
