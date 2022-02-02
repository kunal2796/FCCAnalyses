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

               #MC particles
               .Define("MC_px",  "MCParticle::get_px(Particle)")
               .Define("MC_py",  "MCParticle::get_py(Particle)")
               .Define("MC_pz",  "MCParticle::get_pz(Particle)")
               .Define("MC_e",   "MCParticle::get_e(Particle)")
               .Define("MC_pdg", "MCParticle::get_pdg(Particle)")
               .Define("MC_status","MCParticle::get_genStatus(Particle)")
               
               #reco particles
               .Define("RP_px",  "ReconstructedParticle::get_px(ReconstructedParticles)")
               .Define("RP_py",  "ReconstructedParticle::get_py(ReconstructedParticles)")
               .Define("RP_pz",  "ReconstructedParticle::get_pz(ReconstructedParticles)")
               .Define("RP_e",   "ReconstructedParticle::get_e(ReconstructedParticles)")

               #EE-KT ALGORITHM [status 23]
               #build psedo-jets with the Reconstructed particles
               .Define("pseudo_jets", "JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")
               #add ghosts to pseudo-jets
               #currently selecting status 23 [plan to make it more general soon]
               .Define("pseudo_jets_gm23", "JetClusteringUtils::addGhosts_pseudoJets(pseudo_jets, Particle)")

               #run jet clustering with all reco particles. ee_kt_algorithm, exclusive clustering, exactly 2 jets, E-scheme
               .Define("FCCAnalysesJets_ee_kt23", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pseudo_jets_gm23)")

               #get the jets out of the structure
               .Define("jets_ee_kt23", "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_kt23)")

               #get the jet constituents out of the structure
               .Define("jetconstituents_ee_kt23", "JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_kt23)")

               #get some jet variables
               .Define("jets_ee_kt_e23",  "JetClusteringUtils::get_e(jets_ee_kt23)")
               .Define("jets_ee_kt_px23", "JetClusteringUtils::get_px(jets_ee_kt23)")
               .Define("jets_ee_kt_py23", "JetClusteringUtils::get_py(jets_ee_kt23)")
               .Define("jets_ee_kt_pz23", "JetClusteringUtils::get_pz(jets_ee_kt23)")

               #get jet flavour
               #last argument is the pseudojet vector before adding ghosts [selecting status 23]
               .Define("jets_ee_kt_flavour23", "JetTaggingUtils::get_flavour_gm_auto(jets_ee_kt23, jetconstituents_ee_kt23, Particle, pseudo_jets)")
               #get jet flavour (cone jet definition)
               .Define("jets_ee_kt_flavour_cone","JetTaggingUtils::get_flavour_qqbar(jets_ee_kt23, Particle)")


               #EE-KT ALGORITHM [status 71-79]
               #add ghosts to pseudo-jets
               #currently selecting status 71-79
               .Define("pseudo_jets_gm7x", "JetClusteringUtils::addGhosts7x_pseudoJets(pseudo_jets, Particle)")

               #run jet clustering with all reco particles. ee_kt_algorithm, exclusive clustering, exactly 2 jets, E-scheme
               .Define("FCCAnalysesJets_ee_kt7x", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pseudo_jets_gm7x)")

               #get the jets out of the structure
               .Define("jets_ee_kt7x", "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_kt7x)")

               #get the jet constituents out of the structure
               .Define("jetconstituents_ee_kt7x", "JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_kt7x)")

               #get some jet variables
               .Define("jets_ee_kt_e7x",  "JetClusteringUtils::get_e(jets_ee_kt7x)")
               .Define("jets_ee_kt_px7x", "JetClusteringUtils::get_px(jets_ee_kt7x)")
               .Define("jets_ee_kt_py7x", "JetClusteringUtils::get_py(jets_ee_kt7x)")
               .Define("jets_ee_kt_pz7x", "JetClusteringUtils::get_pz(jets_ee_kt7x)")

               #get jet flavour
               #last argument is the pseudojet vector before adding ghosts [selecting status 7x]
               .Define("jets_ee_kt_flavour7x", "JetTaggingUtils::get_flavour_gm7x_auto(jets_ee_kt7x, jetconstituents_ee_kt7x, Particle, pseudo_jets)")
               
        )

        


        # select branches for output file
        branchList = ROOT.vector('string')()
        for branchName in [
                "MC_px",
                "MC_py",
                "MC_pz",
                "MC_e",
                "MC_pdg",
                "MC_status",

                "RP_px",
                "RP_py",
                "RP_pz",
                "RP_e",

                "jets_ee_kt_e23",
                "jets_ee_kt_px23",
                "jets_ee_kt_py23",
                "jets_ee_kt_pz23",
                "jets_ee_kt_flavour23",
                "jetconstituents_ee_kt23",

                "jets_ee_kt_e7x",
                "jets_ee_kt_px7x",
                "jets_ee_kt_py7x",
                "jets_ee_kt_pz7x",
                "jets_ee_kt_flavour7x",
                "jetconstituents_ee_kt7x",

                "jets_ee_kt_flavour_cone",
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
