import sys
import ROOT

print ("Load cxx analyzers ... ",)
ROOT.gSystem.Load("libedm4hep")
ROOT.gSystem.Load("libpodio")
ROOT.gSystem.Load("libFCCAnalyses")
ROOT.gSystem.Load("libFCCAnalysesFlavour")

ROOT.gErrorIgnoreLevel = ROOT.kFatal
_edm  = ROOT.edm4hep.ReconstructedParticleData()
_pod  = ROOT.podio.ObjectID()
_fcc  = ROOT.dummyLoader
_bs  = ROOT.dummyLoaderFlavour




print ('edm4hep  ',_edm)
print ('podio    ',_pod)
print ('fccana   ',_fcc)

#
#	This is used to process a file in which the Bs and the Bsbar are forced
#	to decay into Jpsi ( -> mu mu) + Phi ( -> K K )
#	We reconstruct the secondary vertex from the 2 muon and 2 kaon tracks.
#       The example also shows how to retrieve the MC and reco'ed Bs legs,
#       as well as the MC Bs, JP]psi and Phis, with their kinematics.
#
#       Example file: 
#       /eos/experiment/fcc/ee/examples/lowerTriangle/p8_ecm91GeV_Zbb_EvtGen_Bs2JpsiPhi_IDEAtrkCov.root
# 	Note: these events were generated at (0,0,0), i.e.no smearing of the
#	primary vertex.
#
#Filter=""

class analysis():

    #__________________________________________________________
    def __init__(self, inputlist, outname, ncpu):
        self.outname = outname
        if ".root" not in outname:
            self.outname+=".root"

        #ROOT.ROOT.EnableImplicitMT(ncpu)
        ROOT.ROOT.DisableImplicitMT()

        self.df = ROOT.RDataFrame("events", inputlist)
        print (" done")
    #__________________________________________________________
    def run(self):
        df2 = (self.df.Range(1000)	# to test over 1000 events only
#        df2 = (self.df

               .Alias("Particle1", "Particle#1.index")
               .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
               .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")


               # MC event primary vertex
               .Define("MC_PrimaryVertex",  "MCParticle::get_EventPrimaryVertex(21)( Particle )" )

               # the recoParticles corresponding  to the tracks that are primaries, according to MC-matching :
               .Define("MC_PrimaryTracks_RP",  "VertexingUtils::SelPrimaryTracks(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle, MC_PrimaryVertex)" )
               # and the corresponding tracks :
               .Define("MC_PrimaryTracks",  "ReconstructedParticle2Track::getRP2TRK( MC_PrimaryTracks_RP, EFlowTrack_1)" )

               # number of tracks in the event
               .Define("ntracks","ReconstructedParticle2Track::getTK_n(EFlowTrack_1)")

               # Retrieve the decay vertex of all MC particles
               #.Define("MC_DecayVertices",  "MCParticle::get_endPoint( Particle, Particle1)" )

               # --- now, determime the primary (and secondary) tracks without using the MC-matching: (from Emmanuel Perez, see https://github.com/HEP-FCC/FCCAnalyses/commit/56a18139761666ebaba134f0243b4aed3d0a6bff)

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

#               # which of the tracks are primary according to the reco algprithm
               .Define("IsPrimary_based_on_reco",  "VertexFitterSimple::IsPrimary_forTracks( EFlowTrack_1, RecoedPrimaryTracks  )")
#               # which of the tracks are primary according to the MC-matching
               .Define("IsPrimary_based_on_MC",  "VertexFitterSimple::IsPrimary_forTracks( EFlowTrack_1,  MC_PrimaryTracks )")
#
#               # jet clustering (ee-kt) before reconstructing SVs in event
               .Define("RP_px",  "ReconstructedParticle::get_px(ReconstructedParticles)")
               .Define("RP_py",  "ReconstructedParticle::get_py(ReconstructedParticles)")
               .Define("RP_pz",  "ReconstructedParticle::get_pz(ReconstructedParticles)")
               .Define("RP_e",   "ReconstructedParticle::get_e(ReconstructedParticles)")
#               #build psedo-jets with the Reconstructed final particles
               .Define("pseudo_jets", "JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")
#               #run jet clustering with all reco particles. ee_kt_algorithm, exclusive clustering, exactly 2 jets, E-scheme
               .Define("FCCAnalysesJets_ee_kt", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pseudo_jets)")
#               #get the jets out of the structure
               .Define("jets_ee_kt", "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_kt)")
#               #get the jet constituents out of the structure
               .Define("jetconstituents_ee_kt", "JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_kt)")

               ### finding SVs in the event ###
               # Event level
               .Define("SV", "VertexFitterSimple::get_SV_event(ReconstructedParticles, EFlowTrack_1, PrimaryVertexObject, IsPrimary_based_on_reco)") # first interface
               #.Define("SV", "VertexFitterSimple::get_SV_event(ReconstructedParticles, EFlowTrack_1, SecondaryTracks, PrimaryVertexObject)")        # second interface
               
               # From jets
               #.Define("SV", "VertexFitterSimple::get_SV_jets(ReconstructedParticles, EFlowTrack_1, PrimaryVertexObject, IsPrimary_based_on_reco, jets_ee_kt, jetconstituents_ee_kt)")               

               # SV properties
               .Define("SV_position", "VertexingUtils::get_position_SV( SV )")
               .Define("ntracks_SV", "VertexingUtils::get_VertexNtrk(SV.vtx)")
               .Define("n_SV", "VertexingUtils::get_n_SV(SV)")
               .Define("SV_mass", "myUtils::get_Vertex_mass( SV.vtx, ReconstructedParticles )")
               .Define("SV_mass_twoPions", "VertexingUtils::get_invM_pairs(SV.vtx)")
               .Define("SV_mass_allPions", "VertexingUtils::get_invM(SV.vtx)")
               .Define("d2PV", "myUtils::get_Vertex_d2PV(VertexFitterSimple::get_all_vertices(PrimaryVertexObject, SV), 0)")
               .Define("d2PV_min", "myUtils::get_dPV2DV_min(d2PV)")
               .Define("d2PV_max", "myUtils::get_dPV2DV_max(d2PV)")
               .Define("d2PV_ave", "myUtils::get_dPV2DV_ave(d2PV)")               

#               .Filter(Filter)

        )


        # select branches for output file
        branchList = ROOT.vector('string')()
        for branchName in [
                "MC_PrimaryVertex",
                "PrimaryVertex",
#                "n_RecoedPrimaryTracks",

                  #  primary vertex and primary tracks w/o any MC-matching :
                  "IsPrimary_based_on_reco",
#                  "PrimaryVertex",

                # vector of SV vertex objects
                "SV",
                "SV_position",
                "SV_mass",
                "SV_mass_twoPions",
                "SV_mass_allPions",
                "n_SV",
                "ntracks_SV",
                "d2PV",
                "d2PV_min",
                "d2PV_max",
                "d2PV_ave",
                "ntracks",
                ]:
            branchList.push_back(branchName)
        df2.Snapshot("events", self.outname, branchList)



if __name__ == "__main__":

    if len(sys.argv)==1:
        print ("usage:")
        print ("python ",sys.argv[0]," file.root")
        sys.exit(3)
    infile = sys.argv[1]
    #outDir = 'FCCee/'+sys.argv[0].split('/')[1]+'/'
    outDir = './'
    import os
    os.system("mkdir -p {}".format(outDir))
    outfile = outDir+infile.split('/')[-1]
    ncpus = 1
    analysis = analysis(infile, outfile, ncpus)
    analysis.run()

    tf = ROOT.TFile(infile)
    entries = tf.events.GetEntries()
    p = ROOT.TParameter(int)( "eventsProcessed", entries)
    outf=ROOT.TFile(outfile,"UPDATE")
    p.Write()
