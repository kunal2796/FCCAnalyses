import sys
import ROOT

print ("Load cxx analyzers ... ",)
ROOT.gSystem.Load("libedm4hep")
ROOT.gSystem.Load("libpodio")
ROOT.gSystem.Load("libFCCAnalyses")
ROOT.gErrorIgnoreLevel = ROOT.kFatal
_edm  = ROOT.edm4hep.ReconstructedParticleData()
_pod  = ROOT.podio.ObjectID()
_fcc  = ROOT.dummyLoader

print ('edm4hep  ',_edm)
print ('podio    ',_pod)
print ('fccana   ',_fcc)


#
# Example file : 
#    /eos/experiment/fcc/ee/examples/lowerTriangle/p8_ecm91GeV_Zuds_IDEAtrkCov.root
#    ( these events were generated at (0,0,0) i.e. no vertex smearing
#
# Example file from spring2021, with vertex smearing :
#    /eos/experiment/fcc/ee/generation/DelphesEvents/spring2021/IDEA/p8_ee_Zuds_ecm91/events_125841058.root


class analysis():

    #__________________________________________________________
    def __init__(self, inputlist, outname, ncpu):
        self.outname = outname
        if ".root" not in outname:
            self.outname+=".root"

        #ROOT.ROOT.EnableImplicitMT(ncpu)

        self.df = ROOT.RDataFrame("events", inputlist)
        print (" done")
    #__________________________________________________________
    def run(self):
        df2 = (self.df.Range(0,1000)
        #df2 = (self.df

               # MC event primary vertex
               .Define("MC_PrimaryVertex",  "MCParticle::get_EventPrimaryVertex(21)( Particle )" )

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
               .Define("SV_jet", "VertexFitterSimple::get_SV_jets(ReconstructedParticles, EFlowTrack_1, PrimaryVertexObject, IsPrimary_based_on_reco, jets_ee_kt, jetconstituents_ee_kt)")

               # finding SVs in the event
               .Define("SV_evt1", "VertexFitterSimple::get_SV_event(ReconstructedParticles, EFlowTrack_1, PrimaryVertexObject, IsPrimary_based_on_reco")
               .Define("SV_evt2", "VertexFitterSimple::get_SV_event(SecondaryTracks, PrimaryVertexObject")

        )


        # select branches for output file
        branchList = ROOT.vector('string')()
        for branchName in [

		#  primary vertex and primary tracks w/o any MC-matching :
		"IsPrimary_based_on_reco",
		"PrimaryVertex",

                # vector of SV vertex objects
                "SV_jet",
                "SV_evt1",
                "SV_evt2",

                ]:
            branchList.push_back(branchName)
        df2.Snapshot("events", self.outname, branchList)



if __name__ == "__main__":

    if len(sys.argv)==1:
        print ("usage:")
        print ("python ",sys.argv[0]," file.root")
        sys.exit(3)
    infile = sys.argv[1]
    outDir = 'FCCee/'+sys.argv[0].split('/')[1]+'/'
    import os
    os.system("mkdir -p {}".format(outDir))
    outfile = outDir+infile.split('/')[-1]
    ncpus = 0
    analysis = analysis(infile, outfile, ncpus)
    analysis.run()

    tf = ROOT.TFile(infile)
    entries = tf.events.GetEntries()
    p = ROOT.TParameter(int)( "eventsProcessed", entries)
    outf=ROOT.TFile(outfile,"UPDATE")
    p.Write()