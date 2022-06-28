#Mandatory: List of processes
processList = {
    'p8_ee_Zuds_ecm91':{'fraction':0.0001}, #Run 0.01% statistics in one output file named <outputDir>/p8_ee_Zuds_ecm91.root
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/spring2021/IDEA/"

#Optional: output directory, default is local running directory
outputDir   = "outputs/FCCee/KG"

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
            .Alias("Particle1", "Particle#1.index")
            
            # get all the MC particles to check for Ks                                                                                                                    
            .Define("MC_pdg", "FCCAnalyses::MCParticle::get_pdg(Particle)")
            # get momenta & mass of all particles                                                                                                                         
            .Define("MC_p4", "FCCAnalyses::MCParticle::get_tlv(Particle)")
            .Define("MC_mass", "FCCAnalyses::MCParticle::get_mass(Particle)")
            
            # Ks -> pi+pi-                                                                                                                                                
            .Define("K0spipi_indices", "FCCAnalyses::MCParticle::get_indices_ExclusiveDecay(310, {211, -211}, true, true) (Particle, Particle1)")
            # Lambda0 -> p+pi-                                                                                                                                            
            .Define("Lambda0ppi_indices", "FCCAnalyses::MCParticle::get_indices_ExclusiveDecay(3122, {2212, -211}, true, true) (Particle, Particle1)")
            
            # determime the primary (and secondary) tracks without using the MC-matching:                                                                                 
            
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
            
            # find V0s                                                                                                                                                    
            .Define("V0", "VertexFinderLCFIPlus::get_V0s(SecondaryTracks, PrimaryVertexObject)")
            # get pdg vector out                                                                                                                                          
            .Define("V0_pdg", "VertexingUtils::get_pdg_V0(V0)")
            # get invariant mass vector out                                                                                                                               
            .Define("V0_invM", "VertexingUtils::get_invM_V0(V0)")
            # get the position                                                                                                                                            
            .Define("V0_pos", "VertexingUtils::get_position_SV(V0)")
            # get the chi2
            .Define("V0_chi2", "VertexingUtils::get_chi2_SV(V0)")
            # get the momenta                                                                                                                                            
            .Define("V0_p", "VertexingUtils::get_p_SV(V0)")

        )
        return df2

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
            # MC particles                                                                                                                                               
            "MC_pdg",
            "MC_p4",
            "MC_mass",
            
            # Ks -> pi+pi- & Lambda0->p+pi-                                                                                                                              
            "K0spipi_indices",
            "Lambda0ppi_indices",
            
            #  primary vertex and primary tracks w/o any MC-matching :                                                                                                   
            "IsPrimary_based_on_reco",
            "PrimaryVertex",
            
            # V0 object                                                                                                                                                  
            "V0",
            "V0_pdg",
            "V0_invM",
            "V0_pos",
            "V0_chi2",
            "V0_p"
        ]
        return branchList
