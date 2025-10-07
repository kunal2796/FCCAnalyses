# testFile="root://eospublic.cern.ch//eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/p8_ee_Zbb_ecm91/events_066726720.root"

# COMMENTED OUT FOR ee, mumu, and tautau
# sv_track_index_*

#Mandatory: List of processes  
processList = {
    # 'p8_ee_Zuds_ecm91':{'fraction':0.0001}, #Run 0.01% statistics in one output file named <outputDir>/p8_ee_Zuds_ecm91.root
    # 'p8_ee_Zcc_ecm91':{'fraction':0.0001}, #Run 0.01% statistics in one output file named <outputDir>/p8_ee_Zcc_ecm91.root
    # 'p8_ee_Zbb_ecm91':{'fraction':0.0001, 'output':'p8_ee_Zbb_ecm91_SV_100K'},
    # 'p8_ee_Zbb_ecm91':{'fraction':0.0001, 'output':'test', 'nevents':1000},

    # 'p8_ee_Zee_ecm91'    :{'fraction':0.001, 'output':'test', 'nevents':1000},
    # 'p8_ee_Zee_ecm91'    :{'fraction':0.015, 'chunks':15},
    # 'p8_ee_Zmumu_ecm91'  :{'fraction':0.015, 'chunks':15},
    # 'p8_ee_Ztautau_ecm91':{'fraction':0.015, 'chunks':15},

    # 'p8_ee_Zbb_ecm91':{'fraction':0.023, 'chunks':100},
    # 'p8_ee_Zcc_ecm91':{'fraction':0.02,  'chunks':100},
    # 'p8_ee_Zss_ecm91':{'fraction':0.02,  'chunks':100},
    # 'p8_ee_Zud_ecm91':{'fraction':0.041,  'chunks':200},
    # 'p8_ee_Zbb_noFSR_ecm91':{'fraction':1, 'chunks':20},

    # 'wzp6_ee_mumu_ecm91p188'              :{'fraction':0.01, 'chunks':10},
    # 'wzp6_ee_mumu_ecm91p188_noISR'        :{'fraction':0.01, 'chunks':10},
    # 'wzp6_ee_mumu_ecm91p188_noISR_noFSR'  :{'fraction':0.01, 'chunks':10},

    # 'wzp6_ee_uu_ecm91p2':{'fraction':0.02, 'chunks':20},
    # 'wzp6_ee_dd_ecm91p2':{'fraction':0.02, 'chunks':20},
    # 'wzp6_ee_ss_ecm91p2':{'fraction':0.02, 'chunks':20},
    # 'wzp6_ee_cc_ecm91p2':{'fraction':0.02, 'chunks':20},
    # 'wzp6_ee_bb_ecm91p2':{'fraction':0.02, 'chunks':20},

    # 'p8_ee_Zbb_ecm91':{'fraction':0.0002},
    'p8_ee_Zss_ecm91':{'fraction':0.002,  'chunks':10},
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/winter2023/IDEA/"

#Optional: output directory, default is local running directory
# outputDir   = "/afs/cern.ch/work/k/kgautam/private/FCCAnalyses/outputs/FCCee/Zqq/asymmetries"
# outputDir   = "/eos/user/k/kgautam/FCCAnalyses/outputs/asymmetries"
# outputDir   = "/eos/user/k/kgautam/FCCAnalyses/outputs/asymmetries/withoutSV"
outputDir   = "/eos/user/k/kgautam/FCCAnalyses/outputs/jetalgo_comparison"

#Optional: ncpus, default is 4
# nCPUS       = 8
# nCPUS       = 1

#Optional running on HTCondor, default is False
runBatch    = True

#Optional batch queue name when running on HTCondor, default is workday  
#batchQueue = "microcentury"
# batchQueue = "longlunch"
batchQueue = "workday"
# batchQueue = "espresso"

#Optional computing account when running on HTCondor, default is group_u_FCC.local_gen
compGroup = "group_u_FCC.local_gen"

############### additions for flavour tagger #################
import os
import urllib.request
# ____________________________________________________________
def get_file_path(url, filename):
    if os.path.exists(filename):
        return os.path.abspath(filename)
    else:
        urllib.request.urlretrieve(url, os.path.basename(url))
        return os.path.basename(url)
# ____________________________________________________________
#
## latest particle transformer model, trainied on 9M jets in winter2023 samples
model_name = "fccee_flavtagging_edm4hep_wc_v1"

## model files needed for unit testing in CI
url_model_dir = "https://fccsw.web.cern.ch/fccsw/testsamples/jet_flavour_tagging/winter2023/wc_pt_13_01_2022/"
url_preproc = "{}/{}.json".format(url_model_dir, model_name)
url_model = "{}/{}.onnx".format(url_model_dir, model_name)

## model files locally stored on /eos
#model_dir = "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_13_01_2022/"
model_dir = "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_7classes_12_04_2023/"
local_preproc = "{}/{}.json".format(model_dir, model_name)
local_model = "{}/{}.onnx".format(model_dir, model_name)

## get local file, else download from url
weaver_preproc = get_file_path(url_preproc, local_preproc)
weaver_model = get_file_path(url_model, local_model)

from addons.ONNXRuntime.jetFlavourHelper import JetFlavourHelper
from addons.FastJet.jetClusteringHelper import ExclusiveJetClusteringHelper

jetFlavourHelper = None
jetClusteringHelper = None
##############################################################

#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):

        ##############################################################
        global jetClusteringHelper
        global jetFlavourHelper

        from examples.FCCee.weaver.config import collections, njets

        tag = ""

        ## define jet clustering parameters
        jetClusteringHelper = ExclusiveJetClusteringHelper(collections["PFParticles"], njets, tag)

        ## run jet clustering
        df = jetClusteringHelper.define(df)

        ## define jet flavour tagging parameters

        jetFlavourHelper = JetFlavourHelper(
            collections,
            jetClusteringHelper.jets,
            jetClusteringHelper.constituents,
            tag,
        )

        ## define observables for tagger
        df = jetFlavourHelper.define(df)

        ## tagger inference
        df = jetFlavourHelper.inference(weaver_preproc, weaver_model, df)
        ##############################################################
        
        df2 = (df

               # alias
               .Alias("Particle0", "Particle#0.index") 
               .Alias("Particle1", "Particle#1.index") 
               .Alias("Jet3","Jet#3.index")

               # Reco'd particles
               .Define("RP_px",        "ReconstructedParticle::get_px(ReconstructedParticles)")
               .Define("RP_py",        "ReconstructedParticle::get_py(ReconstructedParticles)")
               .Define("RP_pz",        "ReconstructedParticle::get_pz(ReconstructedParticles)")               
               .Define("RP_e",         "ReconstructedParticle::get_e(ReconstructedParticles)")

               #build psedo-jets with the Reconstructed final particles
               .Define("pseudo_jets", "JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")

               #ANTI-KT ALGORITHM
               #run jet clustering with all MC particles. kt_algorithm, R=0.5, exclusive clustering, exactly 4 jets, E-scheme
               .Define("FCCAnalysesJets_antikt", "JetClustering::clustering_antikt(0.5, 0, 0, 1, 0)(pseudo_jets)")
               #get the jets out of the structure
               .Define("jets_antikt", "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_antikt)")
               #get the jet constituents out of the structure
               # .Define("jetconstituents_antikt", "JetClusteringUtils::get_constituents(FCCAnalysesJets_antikt)")
               #get some jet variables
               .Define("jets_antikt_e",  "JetClusteringUtils::get_e(jets_antikt)")
               .Define("jets_antikt_px", "JetClusteringUtils::get_px(jets_antikt)")
               .Define("jets_antikt_py", "JetClusteringUtils::get_py(jets_antikt)")
               .Define("jets_antikt_pz", "JetClusteringUtils::get_pz(jets_antikt)")

               #EE-GENKT ALGORITHM
               #run jet clustering with all reconstructed particles. ee_genkt_algorithm, R=0.5, inclusive clustering, E-scheme 
               .Define("FCCAnalysesJets_ee_genkt", "JetClustering::clustering_ee_genkt(0.5, 0, 0, 1, 0, 1)(pseudo_jets)")
               #get the jets out of the struct
               .Define("jets_ee_genkt",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_genkt)")
               #get the jets constituents out of the struct
               # .Define("jetconstituents_ee_genkt","JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_genkt)")
               #get some variables
               .Define("jets_ee_genkt_e",         "JetClusteringUtils::get_e(jets_ee_genkt)")
               .Define("jets_ee_genkt_px",        "JetClusteringUtils::get_px(jets_ee_genkt)")
               .Define("jets_ee_genkt_py",        "JetClusteringUtils::get_py(jets_ee_genkt)")
               .Define("jets_ee_genkt_pz",        "JetClusteringUtils::get_pz(jets_ee_genkt)")

               #EE-KT ALGORITHM
               # jet clustering (ee-kt) before reconstructing SVs in event
               #run jet clustering with all reco particles. ee_kt_algorithm, exclusive clustering, exactly 2 jets, E-scheme
               .Define("FCCAnalysesJets_ee_kt", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pseudo_jets)")
               #get the jets out of the structure
               .Define("jets_ee_kt", "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_kt)")
               #get the jet constituents out of the structure
	       # .Define("jetconstituents_ee_kt", "JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_kt)")
               #get some variables
               # .Define("jets_ee_kt_p",          "JetClusteringUtils::get_p(jets_ee_kt)")
               .Define("jets_ee_kt_px",         "JetClusteringUtils::get_px(jets_ee_kt)")
               .Define("jets_ee_kt_py",         "JetClusteringUtils::get_py(jets_ee_kt)")
               .Define("jets_ee_kt_pz",         "JetClusteringUtils::get_pz(jets_ee_kt)")
               # .Define("jets_ee_kt_theta",      "JetClusteringUtils::get_theta(jets_ee_kt)")
               # .Define("jets_ee_kt_phi",        "JetClusteringUtils::get_phi(jets_ee_kt)")
               .Define("jets_ee_kt_e",          "JetClusteringUtils::get_e(jets_ee_kt)")
               # .Define("n_jets_ee_kt",          "jets_ee_kt.size()")
                             
               )
        return df2

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [

            # EE-KT
            "jets_ee_kt_px",
            "jets_ee_kt_py",
            "jets_ee_kt_pz",
            "jets_ee_kt_e",

            # EE-GENKT
            "jets_ee_genkt_px",
            "jets_ee_genkt_py",
            "jets_ee_genkt_pz",
            "jets_ee_genkt_e",

            # ANTI-KT
            "jets_antikt_px",
            "jets_antikt_py",
            "jets_antikt_pz",
            "jets_antikt_e",

                ]

        # ##  outputs jet properties
        # branchList += jetClusteringHelper.outputBranches()

        # ## outputs jet scores and constituent breakdown
        # branchList += jetFlavourHelper.outputBranches()

        return branchList
