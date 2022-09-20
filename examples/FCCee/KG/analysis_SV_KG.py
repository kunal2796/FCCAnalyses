#Mandatory: List of processes
processList = {
    'p8_ee_Zuds_ecm91':{'fraction':0.0001}, #Run 0.01% statistics in one output file named <outputDir>/p8_ee_Zuds_ecm91.root
    #'p8_ee_Zcc_ecm91':{'fraction':0.0001}, #Run 0.01% statistics in one output file named <outputDir>/p8_ee_Zcc_ecm91.root
    #'p8_ee_Zbb_ecm91':{'fraction':0.0001, 'output':'p8_ee_Zbb_ecm91_SV_100K'},

    #'p8_ee_Zuds_ecm91':{'fraction':0.0003, 'chunks':3}, #Run 0.03% statistics in one output file named <outputDir>/p8_ee_Zuds_ecm91.root
    #'p8_ee_Zbb_ecm91':{'fraction':0.0001, 'chunks':1},
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/spring2021/IDEA/"

#Optional: output directory, default is local running directory
outputDir   = "outputs/FCCee/KG/"

#Optional: ncpus, default is 4
nCPUS       = 1

#Optional running on HTCondor, default is False
#runBatch    = True

#Optional batch queue name when running on HTCondor, default is workday
#batchQueue = "microcentury"

#Optional computing account when running on HTCondor, default is group_u_FCC.local_gen
#compGroup = "group_u_FCC.local_gen"

#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df2 = (df
               # MC event primary vertex
               .Define("MC_PrimaryVertex",  "FCCAnalyses::MCParticle::get_EventPrimaryVertex(21)( Particle )" )
               
               # number of tracks
               .Define("ntracks","ReconstructedParticle2Track::getTK_n(EFlowTrack_1)")
               
               
               # Select primary tracks based on the matching to MC
               # This can be used  to select primary tracks when the
               # gen-level primary vertex  is not  (0,0,0)
               .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
               .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
               # the recoParticles corresponding  to the tracks that are primaries, according to MC-matching :
               .Define("MC_PrimaryTracks_RP",  "VertexingUtils::SelPrimaryTracks(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle, MC_PrimaryVertex)" )
               # and the corresponding tracks :
               .Define("MC_PrimaryTracks",  "ReconstructedParticle2Track::getRP2TRK( MC_PrimaryTracks_RP, EFlowTrack_1)" )
            
               .Define("nPrimaryTracks", "ReconstructedParticle::get_n(MC_PrimaryTracks_RP)")
               
               # Reconstruct the vertex from these primary tracks :
               .Define("VertexObject_primaryTracks",  "VertexFitterSimple::VertexFitter ( 1, MC_PrimaryTracks_RP, EFlowTrack_1) ")  
               .Define("Vertex_primaryTracks",   "VertexingUtils::get_VertexData( VertexObject_primaryTracks )")   # primary vertex, in mm
               
               # Idem, but adding the beam-spot constraint to the vertex fitter
               # At the Z peak, the beam-spot size is ( 4.5 mum, 20 nm, 0.3 mm) 
               # The beam-spot dimensions are passed in mum :
               .Define("VertexObject_primaryTracks_BSC",  "VertexFitterSimple::VertexFitter ( 1, MC_PrimaryTracks_RP, EFlowTrack_1, true, 4.5, 20e-3, 300) ")
               .Define("Vertex_primaryTracks_BSC",   "VertexingUtils::get_VertexData( VertexObject_primaryTracks_BSC )")   # primary vertex, in mm
               
               
               # --- now, determime the primary (and secondary) tracks without using the MC-matching:
               
               # First, reconstruct a vertex from all tracks
               .Define("VertexObject_allTracks",  "VertexFitterSimple::VertexFitter_Tk ( 1, EFlowTrack_1, true, 4.5, 20e-3, 300)")
               # Select the tracks that are reconstructed  as primaries
               .Define("RecoedPrimaryTracks",  "VertexFitterSimple::get_PrimaryTracks( VertexObject_allTracks, EFlowTrack_1, true, 4.5, 20e-3, 300, 0., 0., 0., 0)")
               .Define("n_RecoedPrimaryTracks",  "ReconstructedParticle2Track::getTK_n( RecoedPrimaryTracks )")
               # the final primary vertex : 
               .Define("PrimaryVertexObject",   "VertexFitterSimple::VertexFitter_Tk ( 1, RecoedPrimaryTracks, true, 4.5, 20e-3, 300) ")
               .Define("PrimaryVertex",   "VertexingUtils::get_VertexData( PrimaryVertexObject )")
               
               # the secondary tracks
               .Define("SecondaryTracks",   "VertexFitterSimple::get_NonPrimaryTracks( EFlowTrack_1,  RecoedPrimaryTracks )")
               .Define("n_SecondaryTracks",  "ReconstructedParticle2Track::getTK_n( SecondaryTracks )" )
               
               # which of the tracks are primary according to the reco algprithm
               .Define("IsPrimary_based_on_reco",  "VertexFitterSimple::IsPrimary_forTracks( EFlowTrack_1,  RecoedPrimaryTracks )")
               # which of the tracks are primary according to the MC-matching
               .Define("IsPrimary_based_on_MC",  "VertexFitterSimple::IsPrimary_forTracks( EFlowTrack_1,  MC_PrimaryTracks )")
               
               # jet clustering (ee-kt) before reconstructing SVs in event
               .Define("RP_px",  "ReconstructedParticle::get_px(ReconstructedParticles)")
               .Define("RP_py",  "ReconstructedParticle::get_py(ReconstructedParticles)")
               .Define("RP_pz",  "ReconstructedParticle::get_pz(ReconstructedParticles)")
               .Define("RP_e",   "ReconstructedParticle::get_e(ReconstructedParticles)")
               #build psedo-jets with the Reconstructed final particles
               .Define("pseudo_jets", "JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")
               #run jet clustering with all reco particles. ee_kt_algorithm, exclusive clustering, exactly 2 jets, E-scheme
               .Define("FCCAnalysesJets_ee_kt", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pseudo_jets)")
               #get the jets out of the structure
               .Define("jets_ee_kt", "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_kt)")
               #get the jet constituents out of the structure
               .Define("jetconstituents_ee_kt", "JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_kt)")
               
               
               # finding SVs in jets
               .Define("SV_jet", "VertexFinderLCFIPlus::get_SV_jets(ReconstructedParticles, EFlowTrack_1, PrimaryVertexObject, IsPrimary_based_on_reco, jets_ee_kt, jetconstituents_ee_kt)")
               # finding SVs in the event
               #.Define("SV_evt1", "VertexFinderLCFIPlus::get_SV_event(ReconstructedParticles, EFlowTrack_1, PrimaryVertexObject, IsPrimary_based_on_reco)")
               #.Define("SV_evt2", "VertexFinderLCFIPlus::get_SV_event(SecondaryTracks, PrimaryVertexObject)")
               # multiplicity
               #.Define("SV_evt2_n","VertexingUtils::get_n_SV(SV_evt2)")
               #.Define("SV_evt1_n","VertexingUtils::get_n_SV(SV_evt1)")
               .Define("SV_jet_n", "VertexingUtils::get_n_SV(SV_jet)")
               # vertex position
               .Define("SV_jet_position",  "VertexingUtils::get_position_SV( SV_jet )")
               #.Define("SV_evt1_position", "VertexingUtils::get_position_SV( SV_evt1 )")
               #.Define("SV_evt2_position", "VertexingUtils::get_position_SV( SV_evt2 )")
               # vertex chi2
               .Define("SV_jet_chi2",   "VertexingUtils::get_chi2_SV( SV_jet )")
               #.Define("SV_evt1_chi2",  "VertexingUtils::get_chi2_SV( SV_evt1 )")
               #.Define("SV_evt2_chi2",  "VertexingUtils::get_chi2_SV( SV_evt2 )")
               
              )
        return df2

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
                #  primary vertex and primary tracks w/o any MC-matching :
                'IsPrimary_based_on_reco',
                'PrimaryVertex',
            
                # vector of SV vertex objects
                'SV_jet',
                #'SV_evt1',
                #'SV_evt2',
            
                # SV multiplicity
                #'SV_evt1_n',
                #'SV_evt2_n',
                'SV_jet_n',
                
                # SV position
                'SV_jet_position',
                #'SV_evt1_position',
                #'SV_evt2_position',
                
                # SV chi2
                'SV_jet_chi2',
                #'SV_evt1_chi2',
                #'SV_evt2_chi2',
                ]
        return branchList
