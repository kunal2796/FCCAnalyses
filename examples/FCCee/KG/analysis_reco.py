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

class analysis():

    #__________________________________________________________
    def __init__(self, inputlist, outname, ncpu):
        self.outname = outname
        if ".root" not in outname:
            self.outname+=".root"

        ROOT.ROOT.EnableImplicitMT(ncpu)

        self.df = ROOT.RDataFrame("events", inputlist)
        print (" done")
    #__________________________________________________________
    def run(self):
        df2 = (self.df

               .Alias("Particle0", "Particle#0.index")
               .Alias("Jet3", "Jet#3.index") #alias for dealing with # in python

               .Define("RP_px",  "ReconstructedParticle::get_px(ReconstructedParticles)")
               .Define("RP_py",  "ReconstructedParticle::get_py(ReconstructedParticles)")
               .Define("RP_pz",  "ReconstructedParticle::get_pz(ReconstructedParticles)")
               .Define("RP_p",  "ReconstructedParticle::get_p(ReconstructedParticles)")
               .Define("RP_e",  "ReconstructedParticle::get_e(ReconstructedParticles)")
               .Define("RP_theta",  "ReconstructedParticle::get_theta(ReconstructedParticles)")
               .Define("RP_mass",  "ReconstructedParticle::get_mass(ReconstructedParticles)")
               .Define("RP_charge",  "ReconstructedParticle::get_charge(ReconstructedParticles)")

               #track parameters
               .Define("RP_D0", "ReconstructedParticle2Track::geRP2TRKt_D0(ReconstructedParticles, EFlowTrack_1)")               
               .Define("RP_Z0", "ReconstructedParticle2Track::geRP2TRKt_Z0(ReconstructedParticles, EFlowTrack_1)")               
               .Define("RP_phi", "ReconstructedParticle2Track::geRP2TRKt_phi(ReconstructedParticles, EFlowTrack_1)")               
               .Define("RP_omega", "ReconstructedParticle2Track::geRP2TRKt_omega(ReconstructedParticles, EFlowTrack_1)")               
               .Define("RP_tanLambda", "ReconstructedParticle2Track::geRP2TRKt_tanLambda(ReconstructedParticles, EFlowTrack_1)")               
               
               #track parameters
               .Define("RP_D0_sig", "ReconstructedParticle2Track::geRP2TRKt_D0_sig(ReconstructedParticles, EFlowTrack_1)")               
               .Define("RP_Z0_sig", "ReconstructedParticle2Track::geRP2TRKt_Z0_sig(ReconstructedParticles, EFlowTrack_1)")               

               
               #KT ALGORITHM
               #build psedo-jets with the Reconstructed final particles
               .Define("pseudo_jets", "JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")

               #run jet clustering with all MC particles. kt_algorithm, R=0.5, exclusive clustering, exactly 4 jets, E-scheme
               .Define("FCCAnalysesJets_kt", "JetClustering::clustering_kt(0.5, 2, 2, 1, 0)(pseudo_jets)")

               #get the jets out of the structure
               .Define("jets_kt", "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_kt)")

               #get the jet constituents out of the structure
               .Define("jetconstituents_kt", "JetClusteringUtils::get_constituents(FCCAnalysesJets_kt)")

               #get some jet variables
               .Define("jets_kt_e",  "JetClusteringUtils::get_e(jets_kt)")
               .Define("jets_kt_px", "JetClusteringUtils::get_px(jets_kt)")
               .Define("jets_kt_py", "JetClusteringUtils::get_py(jets_kt)")
               .Define("jets_kt_pz", "JetClusteringUtils::get_pz(jets_kt)")

               #get jet flavour
               .Define("jets_kt_flavour", "JetTaggingUtils::get_flavour(jets_kt)")

               #EE-GENKT ALGORITHM
               #run jet clustering with all reconstructed particles. ee_genkt_algorithm, R=0.5, inclusive clustering, E-scheme 
               .Define("FCCAnalysesJets_ee_genkt", "JetClustering::clustering_ee_genkt(0.5, 0, 0, 1, 0, 1)(pseudo_jets)")

               #get the jets out of the struct
               .Define("jets_ee_genkt",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_genkt)")

               #get the jets constituents out of the struct
               .Define("jetconstituents_ee_genkt","JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_genkt)")

               #get some variables
               .Define("jets_ee_genkt_e",         "JetClusteringUtils::get_e(jets_ee_genkt)")
               .Define("jets_ee_genkt_px",        "JetClusteringUtils::get_px(jets_ee_genkt)")
               .Define("jets_ee_genkt_py",        "JetClusteringUtils::get_py(jets_ee_genkt)")
               .Define("jets_ee_genkt_pz",        "JetClusteringUtils::get_pz(jets_ee_genkt)")

               #get jet flavour
               .Define("jets_ee_genkt_flavour", "JetTaggingUtils::get_flavour(jets_ee_genkt)")

               
               #EE-KT ALGORITHM
               #run jet clustering with all MC particles. ee_kt_algorithm, exclusive clustering, exactly 2 jets, E-scheme
               .Define("FCCAnalysesJets_ee_kt", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pseudo_jets)")

               #get the jets out of the structure
               .Define("jets_ee_kt", "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_kt)")

               #get the jet constituents out of the structure
               .Define("jetconstituents_ee_kt", "JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_kt)")

               #get some jet variables
               .Define("jets_ee_kt_e",  "JetClusteringUtils::get_e(jets_ee_kt)")
               .Define("jets_ee_kt_px", "JetClusteringUtils::get_px(jets_ee_kt)")
               .Define("jets_ee_kt_py", "JetClusteringUtils::get_py(jets_ee_kt)")
               .Define("jets_ee_kt_pz", "JetClusteringUtils::get_pz(jets_ee_kt)")

               #get jet flavour
               .Define("jets_ee_kt_flavour", "JetTaggingUtils::get_flavour(jets_ee_kt)")
               
        )

        


        # select branches for output file
        branchList = ROOT.vector('string')()
        for branchName in [
                "RP_px",
                "RP_py",
                "RP_pz",
                "RP_p",
                "RP_e",
                "RP_theta",
                "RP_mass",
                "RP_charge",

                "RP_D0",
                "RP_Z0",
                "RP_phi",
                "RP_omega",
                "RP_tanLambda",

                "RP_D0_sig",
                "RP_Z0_sig",
                
                "jets_kt_e",
                "jets_kt_px",
                "jets_kt_py",
                "jets_kt_pz",
                "jets_kt_flavour",
                "jetconstituents_kt",
                
                "jets_ee_genkt_e",
                "jets_ee_genkt_px",
                "jets_ee_genkt_py",
                "jets_ee_genkt_pz",
                "jets_ee_genkt_flavour",
                "jetconstituents_ee_genkt",
                
                "jets_ee_kt_e",
                "jets_ee_kt_px",
                "jets_ee_kt_py",
                "jets_ee_kt_pz",
                "jets_ee_kt_flavour",
                "jetconstituents_ee_kt",
                ]:
            branchList.push_back(branchName)
        df2.Snapshot("events", self.outname, branchList)

# example call for standalone file
# python examples/FCCee/higgs/mH-recoil/mumu/analysis.py /eos/experiment/fcc/ee/generation/DelphesEvents/fcc_tmp/p8_ee_ZH_ecm240/events_058720051.root
if __name__ == "__main__":

    if len(sys.argv)==1:
        print ("usage:")
        print ("python ",sys.argv[0]," file.root")
        sys.exit(3)
    infile = sys.argv[1]
    outDir = sys.argv[0].replace(sys.argv[0].split('/')[0],'outputs').replace('analysis.py','')+'/'
    import os
    os.system("mkdir -p {}".format(outDir))
    outfile = outDir+infile.split('/')[-1]
    ncpus = 0
    print ('outfile  ',outfile)
    analysis = analysis(infile, outfile, ncpus)
    analysis.run()

    tf = ROOT.TFile(infile)
    entries = tf.events.GetEntries()
    p = ROOT.TParameter(int)( "eventsProcessed", entries)
    outf=ROOT.TFile(outfile,"UPDATE")
    p.Write()
