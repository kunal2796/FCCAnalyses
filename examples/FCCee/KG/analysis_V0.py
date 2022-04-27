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
        df2 = (self.df.Range(0,10000)
        #df2 = (self.df

               .Alias("Particle1", "Particle#1.index")

               # get all the MC particles to check for Ks
               .Define("MC_pdg", "MCParticle::get_pdg(Particle)")
               # get momenta & mass of all particles
               .Define("MC_p4", "MCParticle::get_tlv(Particle)")
               .Define("MC_mass", "MCParticle::get_mass(Particle)")

               # Ks -> pi+pi-
               .Define("K0spipi_indices", "MCParticle::get_indices_ExclusiveDecay(310, {211, -211}, true, true) (Particle, Particle1)")
               # Lambda0 -> p+pi-
               .Define("Lambda0ppi_indices", "MCParticle::get_indices_ExclusiveDecay(3122, {2212, -211}, true, true) (Particle, Particle1)")

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
               .Define("V0", "VertexFitterSimple::get_V0s(SecondaryTracks, PrimaryVertexObject)")
               # get pdg vector out
               .Define("V0_pdg", "VertexingUtils::get_pdg_V0(V0)")
               # get invariant mass vector out
               .Define("V0_invM", "VertexingUtils::get_invM_V0(V0)")
               # get the position
               .Define("V0_pos", "VertexingUtils::get_position_SV(V0)")

        )


        # select branches for output file
        branchList = ROOT.vector('string')()
        for branchName in [

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
