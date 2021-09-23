import timeit
import uproot
import awkward as ak
import os
import h5py as h5
import numpy as np
import sparse
import matplotlib.pyplot as plt
from collections import defaultdict
from coffea import hist, processor
from coffea.nanoevents import NanoEventsFactory, BaseSchema, NanoAODSchema
from coffea.nanoevents.methods import vector

#import sys
#sys.modules['__main__'].__file__ = 'ipython'

ak.behavior.update(vector.behavior)

# https://stackoverflow.com/questions/1987694/how-to-print-the-full-numpy-array-without-truncation/24542498#24542498
def fullprint(*args, **kwargs):
  from pprint import pprint
  import numpy
  opt = numpy.get_printoptions()
  numpy.set_printoptions(threshold=numpy.inf)
  pprint(*args, **kwargs)
  numpy.set_printoptions(**opt)

def len_from_pt(pt):
    #avoid calling as this is not array-at-a-time
    return [len(x) for x in pt] 

def Writeh5(output, name,folder):
    if not os.path.exists(folder):
        os.makedirs(folder)
    with h5.File(os.path.join(folder,'{0}.h5'.format(name)), "w") as fh5:
        dset = fh5.create_dataset("data", data=output['data'].value)
        dset = fh5.create_dataset("data_len", data=output['data_len'].value)
        dset = fh5.create_dataset("coords", data=output['coords'].value)
        dset = fh5.create_dataset("shape", data=output['shape'].value)
        dset = fh5.create_dataset("fill_value", data=output['fill_value'].value)
        #dset = fh5.create_dataset("data", data=output['data'])
        dset = fh5.create_dataset("pid", data=output['pid'].value)
        #dset = fh5.create_dataset("global", data=output['global'])
        #dset = fh5.create_dataset("px", data=output['px'])
        #dset = fh5.create_dataset("py", data=output['py'])
        #dset = fh5.create_dataset("pz", data=output['pz'])
        #dset = fh5.create_dataset("E", data=output['E'])
        #dset = fh5.create_dataset("invMass", data=output['invMass'])


class Channel():
    def __init__(self, name):
        self.name = name
        self.sel=None
        self.jets= None
        self.pfs=None
        self.pfs_pid=None
        self.glob = None
        self.jet = None
        self.pid = None
        

    def apply_jet_sel(self,quark='up'):
        jet_matched_pfs = self.jets[self.pfs.idx]

        #print(gen_electrons.pdgId)
        #Cuts on jets

        if quark=='up':
            pfs_mask = (np.abs(jet_matched_pfs.partonFlavour)==2) & (jet_matched_pfs.hadronFlavour==0)
        elif quark == 'strange':
            pfs_mask = (np.abs(jet_matched_pfs.partonFlavour)==3) & (jet_matched_pfs.hadronFlavour==0)
        elif quark == 'down':
            pfs_mask = (np.abs(jet_matched_pfs.partonFlavour)==1) & (jet_matched_pfs.hadronFlavour==0)
        else:
            pfs_mask = (np.abs(jet_matched_pfs.partonFlavour)==21) & (jet_matched_pfs.hadronFlavour==0)
            
        self.pfs = self.pfs[pfs_mask]


class MyProcessor(processor.ProcessorABC):
    def __init__(self):
        print('ONCE')
        
        #self._accumulator = processor.dict_accumulator({
        #    "sumw": processor.list_accumulator(),
        #})
        #self._accumulator = processor.list_accumulator()
        
        ##dict_accumulator = {}
        ML_accumulator = {
            'data': processor.column_accumulator(np.zeros([0,])),
            'data_len': processor.column_accumulator(np.zeros([0], dtype=int)),
            'coords': processor.column_accumulator(np.zeros([0,4], dtype=int)),
            'shape': processor.column_accumulator(np.zeros([0], dtype=int)),
            'fill_value': processor.column_accumulator(np.zeros([0])),
            'pid': processor.column_accumulator(np.zeros([0])),
        }
        self._accumulator = processor.dict_accumulator( ML_accumulator )
        '''
        ML_accumulator = {
            
            'data': processor.list_accumulator(),
            'px': processor.list_accumulator(),
            'py': processor.list_accumulator(),
            'pz': processor.list_accumulator(),
            'E': processor.list_accumulator(),
            'pid': processor.list_accumulator(),
            'MC_px': processor.list_accumulator(),
            'MC_py': processor.list_accumulator(),
            'MC_pz': processor.list_accumulator(),
            'MC_E': processor.list_accumulator(),
            'MC_E_': processor.list_accumulator(),
            'MC_pdg': processor.list_accumulator(),
            'invMass': processor.list_accumulator(),
            'global': processor.list_accumulator(),
        }
        '''
        ##dict_accumulator['ML'] = processor.dict_accumulator(ML_accumulator)
        ##self._accumulator = processor.dict_accumulator( dict_accumulator )

    @property
    def accumulator(self):
        return self._accumulator



    def process(self, events):
        output = self.accumulator.identity()
        
        jets = ak.zip({
            "t": events.jets_ee_kt_e,
            "x": events.jets_ee_kt_px,
            "y": events.jets_ee_kt_py,
            "z": events.jets_ee_kt_pz,
            "partonFlavour": events.jets_ee_kt_flavour,
            #"r": np.sqrt(events.jets_ee_kt_px**2+events.jets_ee_kt_py**2),
            #"theta": np.arctan2(np.sqrt(events.jets_ee_kt_px**2+events.jets_ee_kt_py**2), events.jets_ee_kt_pz),
        ##}, with_name="PtEtaPhiMCandidate")
        }, with_name="LorentzVector")
        

        GenPartsf = ak.zip({
            "t": events.MC_e_final,
            "x": events.MC_px_final,
            "y": events.MC_py_final,
            "z": events.MC_pz_final,
            "pdg": events.MC_pdg_final,
            #"r": np.sqrt(events.MC_px_final**2+events.MC_py_final**2),
        }, with_name="LorentzVector")
        print('NEUTRINO TESTS')
        print('GenPartsf')
        print(GenPartsf.pdg)
        GenParts = ak.zip({
            "t": events.MC_e,
            "x": events.MC_px,
            "y": events.MC_py,
            "z": events.MC_pz,
            "pdg": events.MC_pdg,
            #"K0spipi_indices": events.K0spipi_indices,
            "MC_status": events.MC_status,
        }, with_name="LorentzVector")
        
        
        ################
        #Index Matching#
        ################

        #Below is a built in bug that will fail to run the code if there are more than 250 MC particles in an event. I did this so that one would not clip the event at the tail and miss stable particles in the process...
        
        #We begin by defining array_size which is the first output of ak.where -> flattened array corresponding to outer index where condition (MC_status==1 is true)
        # e.g. MC_status=[[0,1,2,3],[2,3,4],[3,2,1,1,0]] --ak.where--> ([0, 2, 2],[1, 2, 3]) 
        #Hence, array_size gives the indices of the events where there are stable MC particles

        array_size = np.array(ak.where(ak.fill_none(ak.pad_none((GenParts.MC_status==1), 300, clip=False), 0))[0][:])
        
        #array_content gives the incdices of the stable MC particles within each event

        array_content = np.array(ak.where(ak.fill_none(ak.pad_none((GenParts.MC_status==1), 300, clip=False), 0))[1][:])
        
        #array_split_content gives the unflattened version of array_content, in other words array_content is split such that each index is in a subarray corresponding to the event which it references
        #array_content = [0,1,2,4,1,2,3,0,4,...] -> array_split_content = [[0,1,2,4],[1,2,3],[0,4],...] 

        array_split_content = ak.unflatten(array_content, np.unique(array_size, return_counts=True)[1])
        
        array_split_content_padded = np.array(ak.fill_none(ak.pad_none(array_split_content, 200, clip=False), -1))
        
        #here we define mask_padded as the first element of events.K0spipi_indices, in other words the index of the first K0s in each event

        
        #mask_padded_K0spipi = events.K0spipi_indices[ak.num(events.K0spipi_indices)>0][:,0]
        #mask_padded_pi0gammagamma = events.pi0gammagamma_indices[ak.num(events.K0spipi_indices)>0][:,0]
        
        #Here we define mask_padded_raw, which corresponds to every 2,3 element of events.K0spipi_indices, in other words the daughter pions
        # e.g. mask_padded_raw = [78,79,96,97,-999, ..., -999]

        mask_padded_raw_K0spipi = np.delete(np.array(ak.fill_none(ak.pad_none(events.K0spipi_indices, 30, clip=False), -999)), np.arange(0, 30, 3), axis=1)
        mask_padded_raw_pi0gammagamma = np.delete(np.array(ak.fill_none(ak.pad_none(events.pi0gammagamma_indices, 30, clip=False), -999)), np.arange(0, 30, 3), axis=1)

        #Below we use the mask_padded_raw and equate it to array_split_content to create a truth mask
        #Due to broadcasting issues we can actually only do this for one pion per event at a time, across all events, below is the first

        K0spipi_mask = (array_split_content_padded==np.expand_dims(mask_padded_raw_K0spipi.T[0], axis=1))
        pi0gammagamma_mask = (array_split_content_padded==np.expand_dims(mask_padded_raw_pi0gammagamma.T[0], axis=1))
        
        #We then proceed as above for the rest of the pions (max 20), by adding the resulting truth masks

        for i in range(1,len(mask_padded_raw_K0spipi.T)):
            holder_pi = (array_split_content_padded==np.expand_dims(mask_padded_raw_K0spipi.T[i], axis=1))
            holder_gamma = (array_split_content_padded==np.expand_dims(mask_padded_raw_pi0gammagamma.T[i], axis=1))
            K0spipi_mask = np.add(K0spipi_mask, holder_pi)
            pi0gammagamma_mask = np.add(pi0gammagamma_mask, holder_gamma)

        #Finally, we define new_indices as the positions of non zero elements with respect to the array_split_content, in other words the indices of the daughter pions in the MC_status==1 collection that we are clustering
        print('---------------------------------------------------------------')
        print(events.pi0gammagamma_indices[0])

         
        
        structmask = ak.fill_none(ak.pad_none((np.abs(GenPartsf.pdg)>=0), 200, clip=False), False)
        print(ak.type(GenPartsf))
        GenPartsf['K0spipi_mask'] = ak.Array(K0spipi_mask)[structmask]
        GenPartsf['pi0gammagamma_mask'] = ak.Array(pi0gammagamma_mask)[structmask]
        print(ak.type(GenPartsf))
        print(ak.where(GenPartsf.pi0gammagamma_mask[0]))
        
        
         
        print('////////////////////////////')
        print(events.K0spipi_indices[0])
        print(GenParts.x[events.K0spipi_indices][0,1::3])
        a = np.delete(np.array(ak.fill_none(ak.pad_none(GenParts.x[events.K0spipi_indices], 30, clip=False), -999)), np.arange(0, 30, 3), axis=1)
        print(a)
        #print(GenPartsf.x[new_indices_K0spipi][0])
        print(GenPartsf.x[GenPartsf.K0spipi_mask][0,0::2])
        b = np.array(ak.fill_none(ak.pad_none(GenPartsf.x[GenPartsf.K0spipi_mask], 20, clip=False), -999))
        print(b)
        print(a[0]==b[0])
        print(np.array_equal(a,b))
        print(len(a[0]))
        print(len(b[0]))
        
        
        #RANDOMIZE FIRST OR SECOND JET!!!

        #Here we begin with the selection of Kaon_shorts, by selecting every third element from events.K0spipi_indices

        Kaon_short = ak.copy(GenParts[events.K0spipi_indices[:,::3]])
        neutralPion = ak.copy(GenParts[events.pi0gammagamma_indices[:,::3]])

        #For unstable particles it is required that the full MC collection be referenced, and as such that angles are computed separately. I begin here with the K0s->pipi by defining angular differences with respect to each jet.
        
        shortKaon_theta_diff1 = jets.theta[:,0] - Kaon_short.theta
        shortKaon_theta_diff2 = jets.theta[:,1] - Kaon_short.theta
        
        shortKaon_phi_diff1 = jets.phi[:,0] - Kaon_short.phi
        shortKaon_phi_diff1 = ((shortKaon_phi_diff1+np.pi)%2)-np.pi
        shortKaon_phi_diff2 = jets.phi[:,1] - Kaon_short.phi
        shortKaon_phi_diff2 = ((shortKaon_phi_diff2+np.pi)%2)-np.pi
        
        neutralPion_theta_diff1 = jets.theta[:,0] - neutralPion.theta
        neutralPion_theta_diff2 = jets.theta[:,1] - neutralPion.theta
        
        neutralPion_phi_diff1 = jets.phi[:,0] - neutralPion.phi
        neutralPion_phi_diff1 = ((neutralPion_phi_diff1+np.pi)%2)-np.pi
        neutralPion_phi_diff2 = jets.phi[:,1] - neutralPion.phi
        neutralPion_phi_diff2 = ((neutralPion_phi_diff2+np.pi)%2)-np.pi
        
        #Next the sum of the angle squares as a proxy to the angle between the two vectors, in order to select the Kaon shorts corresponding to the given jet (only the first in our case)
        
        vec_angle1 = (shortKaon_theta_diff1**2+shortKaon_phi_diff1**2)
        vec_angle2 = (shortKaon_theta_diff2**2+shortKaon_phi_diff2**2)

        vec_angle1_neutralPion = (neutralPion_theta_diff1**2+neutralPion_phi_diff1**2)
        vec_angle2_neutralPion = (neutralPion_theta_diff2**2+neutralPion_phi_diff2**2)
        
        #Note that henceforth shortKaon will be treated as an additional category

        shortKaon_theta_diff1 = shortKaon_theta_diff1[vec_angle1<vec_angle2]
        shortKaon_phi_diff1 = shortKaon_phi_diff1[vec_angle1<vec_angle2]
        
        neutralPion_theta_diff1 = neutralPion_theta_diff1[vec_angle1_neutralPion<vec_angle2_neutralPion]
        neutralPion_phi_diff1 = neutralPion_phi_diff1[vec_angle1_neutralPion<vec_angle2_neutralPion]


        #Note that for now the implementation is hacky and we basically pass the jet constituents mask first for one jet and then for the other (this may become problematic if there are multiple jets), but this is fine for now as we focus on the exclusive ee_kt. In the future we might want to think about how we could pass both jets at once (19/08/2021)...
        
        
        #Here we simply defined the quark mask for each case, used to define the pid, refer to FCCAnalyses for the jet flavour algorithm (should we continue using this even in the case of Z->jj?)

        up_mask = (np.abs(jets.partonFlavour[:,0])==2)&(np.abs(jets.partonFlavour[:,1])==2)&((jets.partonFlavour[:,0]*jets.partonFlavour[:,1])<0)
        down_mask = (np.abs(jets.partonFlavour[:,0])==1)&(np.abs(jets.partonFlavour[:,1])==1)&((jets.partonFlavour[:,0]*jets.partonFlavour[:,1])<0)
        strange_mask = (np.abs(jets.partonFlavour[:,0])==3)&(np.abs(jets.partonFlavour[:,1])==3)&((jets.partonFlavour[:,0]*jets.partonFlavour[:,1])<0)
        charm_mask = (np.abs(jets.partonFlavour[:,0])==4)&(np.abs(jets.partonFlavour[:,1])==4)&((jets.partonFlavour[:,0]*jets.partonFlavour[:,1])<0)
        bottom_mask = (np.abs(jets.partonFlavour[:,0])==5)&(np.abs(jets.partonFlavour[:,1])==5)&((jets.partonFlavour[:,0]*jets.partonFlavour[:,1])<0)

        pid = 2*up_mask+down_mask+3*strange_mask+4*charm_mask+5*bottom_mask
        
        
        #Here firstjetGenPartsf is defined as the collection of stable MC particles corresponding to the first jet

        firstjetGenPartsf = ak.copy(GenPartsf[events.jetconstituents_ee_kt[:,0]])
        

        #more attempts to free up memory...
        secondjetGenPartsf = ak.copy(GenPartsf[events.jetconstituents_ee_kt[:,1]])
        del GenPartsf
        
        print('INDEX CHECK')
        print(ak.where(K0spipi_mask[0]))
        print((events.jetconstituents_ee_kt[0,0]))
        print((events.jetconstituents_ee_kt[0,1]))

        #Here the different categories are defined, corresponding to the different (optimistically) distinguishable particle types 

        neutralKaon_mask = (np.abs(firstjetGenPartsf.pdg)==130)|(np.abs(firstjetGenPartsf.pdg)==310)
        chargedKaon_mask = (np.abs(firstjetGenPartsf.pdg)==321)

        #Worthwhile REMOVING the neutral Pion category (from here)!

        neutralPion_mask = (np.abs(firstjetGenPartsf.pdg)==111)
        chargedPion_mask = (np.abs(firstjetGenPartsf.pdg)==211)&(~firstjetGenPartsf.K0spipi_mask)
        electron_mask = (np.abs(firstjetGenPartsf.pdg)==11)
        muon_mask = (np.abs(firstjetGenPartsf.pdg)==13)
        photon_mask = (np.abs(firstjetGenPartsf.pdg)==22)&(~firstjetGenPartsf.pi0gammagamma_mask)
        protonNeutron_mask = (np.abs(firstjetGenPartsf.pdg)==2212)|(np.abs(firstjetGenPartsf.pdg)==2112)
        
        '''
        kaon_mask = ak.copy(np.abs(firstjetGenPartsf.pdg)==321)
        kaon_parts = firstjetGenPartsf[kaon_mask]
        kaons = (jets[ak.num(kaon_parts.t)>0])[:,0]
        kaons = jets[:,0]
        #kaon_parts = kaon_parts[ak.num(kaon_parts.t)>0]
        '''


        #Here the angular differences for each stable MC particle and the jet is computed

        theta_diff = jets.theta[:,0] - firstjetGenPartsf.theta
        phi_diff = jets.phi[:,0] - firstjetGenPartsf.phi
        phi_diff = ((phi_diff+np.pi)%2)-np.pi

        '''
        #histogramming in a single event seems simple enough, but now I must proceed and apply it to all while avoiding loops cause speed. Helpful suggestion found here: https://stackoverflow.com/questions/18851471/numpy-histogram-on-multi-dimensional-array 
        #side note-> np is fine as the array is jagged overall but not per entry, could use coffea if I want to also (as coffea hist has values() to get values)
        '''
        
        #Here nphisto is defined to hold all histograms (9 classes, 29x29 resolution), for now -0.5->0.5 rad for phi, theta

        nphisto = np.zeros([len(firstjetGenPartsf.phi),9,29,29])
        bins = np.linspace(-0.5,0.5,30)
        
        neutralKaon_phi_diff = phi_diff[neutralKaon_mask]
        chargedKaon_phi_diff = phi_diff[chargedKaon_mask]
        neutralPion_phi_diff = phi_diff[neutralPion_mask]
        chargedPion_phi_diff = phi_diff[chargedPion_mask]
        electron_phi_diff = phi_diff[electron_mask]
        muon_phi_diff = phi_diff[muon_mask]
        photon_phi_diff = phi_diff[photon_mask]
        protonNeutron_phi_diff = phi_diff[protonNeutron_mask]
        
        neutralKaon_theta_diff = theta_diff[neutralKaon_mask]
        chargedKaon_theta_diff = theta_diff[chargedKaon_mask]
        neutralPion_theta_diff = theta_diff[neutralPion_mask]
        chargedPion_theta_diff = theta_diff[chargedPion_mask]
        electron_theta_diff = theta_diff[electron_mask]
        muon_theta_diff = theta_diff[muon_mask]
        photon_theta_diff = theta_diff[photon_mask]
        protonNeutron_theta_diff = theta_diff[protonNeutron_mask]
        
        neutralKaon_weights = ((firstjetGenPartsf.p)/(jets.p[:,0]))[neutralKaon_mask]
        chargedKaon_weights = ((firstjetGenPartsf.p)/(jets.p[:,0]))[chargedKaon_mask]
        neutralPion_weights = ((firstjetGenPartsf.p)/(jets.p[:,0]))[neutralPion_mask]
        chargedPion_weights = ((firstjetGenPartsf.p)/(jets.p[:,0]))[chargedPion_mask]
        electron_weights = ((firstjetGenPartsf.p)/(jets.p[:,0]))[electron_mask]
        muon_weights = ((firstjetGenPartsf.p)/(jets.p[:,0]))[muon_mask]
        photon_weights = ((firstjetGenPartsf.p)/(jets.p[:,0]))[photon_mask]
        protonNeutron_weights = ((firstjetGenPartsf.p)/(jets.p[:,0]))[protonNeutron_mask]

        shortKaon_weights = ((Kaon_short.p[vec_angle1<vec_angle2])/(jets.p[:,0]))
        
        neutralPion_weights = ((neutralPion.p[vec_angle1_neutralPion<vec_angle2_neutralPion])/(jets.p[:,0]))
        
        #Can shortened from 30, esp. for things like muons where we will never have 30 muons per jet, can also set clip to False
        neutralKaon_phi_diff = np.asarray(ak.fill_none(ak.pad_none(neutralKaon_phi_diff, 70, clip=True), 999))
        chargedKaon_phi_diff = np.asarray(ak.fill_none(ak.pad_none(chargedKaon_phi_diff, 70, clip=True), 999))
        neutralPion_phi_diff = np.asarray(ak.fill_none(ak.pad_none(neutralPion_phi_diff, 70, clip=True), 999))
        chargedPion_phi_diff = np.asarray(ak.fill_none(ak.pad_none(chargedPion_phi_diff, 70, clip=True), 999))
        electron_phi_diff = np.asarray(ak.fill_none(ak.pad_none(electron_phi_diff, 70, clip=True), 999))
        muon_phi_diff = np.asarray(ak.fill_none(ak.pad_none(muon_phi_diff, 70, clip=True), 999))
        photon_phi_diff = np.asarray(ak.fill_none(ak.pad_none(photon_phi_diff, 70, clip=True), 999))
        protonNeutron_phi_diff = np.asarray(ak.fill_none(ak.pad_none(protonNeutron_phi_diff, 70, clip=True), 999))
        
        shortKaon_phi_diff = np.asarray(ak.fill_none(ak.pad_none(shortKaon_phi_diff1, 70, clip=True), 999))
        
        neutralPion_phi_diff = np.asarray(ak.fill_none(ak.pad_none(neutralPion_phi_diff1, 70, clip=True), 999))
        #print('SHORT KAON PHI DIFF')
        #print(shortKaon_phi_diff)
        
        neutralKaon_theta_diff = np.asarray(ak.fill_none(ak.pad_none(neutralKaon_theta_diff, 70, clip=True), 999))
        chargedKaon_theta_diff = np.asarray(ak.fill_none(ak.pad_none(chargedKaon_theta_diff, 70, clip=True), 999))
        neutralPion_theta_diff = np.asarray(ak.fill_none(ak.pad_none(neutralPion_theta_diff, 70, clip=True), 999))
        chargedPion_theta_diff = np.asarray(ak.fill_none(ak.pad_none(chargedPion_theta_diff, 70, clip=True), 999))
        electron_theta_diff = np.asarray(ak.fill_none(ak.pad_none(electron_theta_diff, 70, clip=True), 999))
        muon_theta_diff = np.asarray(ak.fill_none(ak.pad_none(muon_theta_diff, 70, clip=True), 999))
        photon_theta_diff = np.asarray(ak.fill_none(ak.pad_none(photon_theta_diff, 70, clip=True), 999))
        protonNeutron_theta_diff = np.asarray(ak.fill_none(ak.pad_none(protonNeutron_theta_diff, 70, clip=True), 999))
        
        shortKaon_theta_diff = np.asarray(ak.fill_none(ak.pad_none(shortKaon_theta_diff1, 70, clip=True), 999))
        
        neutralPion_theta_diff = np.asarray(ak.fill_none(ak.pad_none(neutralPion_theta_diff1, 70, clip=True), 999))
        
        #Write this padding into a function, and a subsequent loop from a list, or maybe channels? Turned out longer and uglier than I expected...
        neutralKaon_weights = np.asarray(ak.fill_none(ak.pad_none(neutralKaon_weights, 70, clip=True), 0))
        chargedKaon_weights = np.asarray(ak.fill_none(ak.pad_none(chargedKaon_weights, 70, clip=True), 0))
        neutralPion_weights = np.asarray(ak.fill_none(ak.pad_none(neutralPion_weights, 70, clip=True), 0))
        chargedPion_weights = np.asarray(ak.fill_none(ak.pad_none(chargedPion_weights, 70, clip=True), 0))
        electron_weights = np.asarray(ak.fill_none(ak.pad_none(electron_weights, 70, clip=True), 0))
        muon_weights = np.asarray(ak.fill_none(ak.pad_none(muon_weights, 70, clip=True), 0))
        photon_weights = np.asarray(ak.fill_none(ak.pad_none(photon_weights, 70, clip=True), 0))
        protonNeutron_weights = np.asarray(ak.fill_none(ak.pad_none(protonNeutron_weights, 70, clip=True), 0))
        
        shortKaon_weights = np.asarray(ak.fill_none(ak.pad_none(shortKaon_weights, 70, clip=True), 0))
        
        neutralPion_weights = np.asarray(ak.fill_none(ak.pad_none(neutralPion_weights, 70, clip=True), 0))
        

        
        # https://iscinumpy.gitlab.io/post/histogram-speeds-in-python/ -> could follow one of these approaches or simply histogram in root entirely... The upside to histogramming in root is that we will have to save generic h5 files which we could feed to any network...
        
        start_time = timeit.default_timer()
        for i in range(len(firstjetGenPartsf.phi)):
            nphisto[i,0] = np.histogram2d(neutralKaon_phi_diff[i], neutralKaon_theta_diff[i], bins=bins, weights=neutralKaon_weights[i])[0]
            nphisto[i,1] = np.histogram2d(chargedKaon_phi_diff[i], chargedKaon_theta_diff[i], bins=bins, weights=chargedKaon_weights[i])[0]
            nphisto[i,2] = np.histogram2d(neutralPion_phi_diff[i], neutralPion_theta_diff[i], bins=bins, weights=neutralPion_weights[i])[0]
            nphisto[i,3] = np.histogram2d(chargedPion_phi_diff[i], chargedPion_theta_diff[i], bins=bins, weights=chargedPion_weights[i])[0]
            nphisto[i,4] = np.histogram2d(electron_phi_diff[i], electron_theta_diff[i], bins=bins, weights=electron_weights[i])[0]
            nphisto[i,5] = np.histogram2d(muon_phi_diff[i], muon_theta_diff[i], bins=bins, weights=muon_weights[i])[0]
            nphisto[i,6] = np.histogram2d(photon_phi_diff[i], photon_theta_diff[i], bins=bins, weights=photon_weights[i])[0]
            nphisto[i,7] = np.histogram2d(protonNeutron_phi_diff[i], protonNeutron_theta_diff[i], bins=bins, weights=protonNeutron_weights[i])[0]
            nphisto[i,8] = np.histogram2d(shortKaon_phi_diff[i], shortKaon_theta_diff[i], bins=bins, weights=shortKaon_weights[i])[0]
        elapsed = timeit.default_timer() - start_time
        print('ELAPSED TIME -> '+str(elapsed))
        
        nphisto = nphisto[pid!=0]
        pid = pid[pid!=0]

        nphisto_sparse = sparse.COO(nphisto)
        output['data']+= processor.column_accumulator(np.asarray(nphisto_sparse.data))#PFCands.pt[0].tolist()
        output['data_len']+= processor.column_accumulator(np.expand_dims(np.asarray(len(nphisto_sparse.data), dtype=int), axis=0))
        output['coords']+= processor.column_accumulator(np.asarray(nphisto_sparse.coords, dtype=int).T)#PFCands.pt[0].tolist()
        output['shape']+= processor.column_accumulator(np.asarray(nphisto_sparse.shape, dtype=int))#PFCands.pt[0].tolist()
        output['fill_value']+= processor.column_accumulator(np.expand_dims(np.asarray(nphisto_sparse.fill_value), axis=0))#PFCands.pt[0].tolist()
        output['pid']+= processor.column_accumulator(np.asarray(pid))

        ##output['sumw']+= jets.partonFlavour.tolist()#PFCands.pt[0].tolist()
        '''
        output['ML']['pid']+= jets.partonFlavour.tolist()#PFCands.pt[0].tolist()
        output['ML']['px']+= jets.x.tolist()#PFCands.pt[0].tolist()
        output['ML']['py']+= jets.y.tolist()#PFCands.pt[0].tolist()
        output['ML']['pz']+= jets.z.tolist()#PFCands.pt[0].tolist()
        output['ML']['E']+= jets.t.tolist()#PFCands.pt[0].tolist()
        output['ML']['invMass']+= jets.mass.tolist()#PFCands.pt[0].tolist()
        output['ML']['MC_pdg_jet']+= GenPartsf.pdg.tolist()#PFCands.pt[0].tolist()
        output['ML']['MC_px_jet']+= GenPartsf.x.tolist()#PFCands.pt[0].tolist()
        output['ML']['MC_py_jet']+= GenPartsf.y.tolist()#PFCands.pt[0].tolist()
        output['ML']['MC_pz_jet']+= GenPartsf.z.tolist()#PFCands.pt[0].tolist()
        #output['ML']['MC_E']+= GenPartsf.t[events.jetconstituents_kt]#PFCands.pt[0].tolist()
        #output['ML']['MC_E_']+= GenPartsf[events.jetconstituents_kt].t.tolist()#PFCands.pt[0].tolist()
        '''
        return output

    def postprocess(self, accumulator):
        return accumulator

#p = MyProcessor()
filename = "rootFiles/p8_ee_Zuds_ecm91.root"
file = uproot.open(filename)
'''
events = NanoEventsFactory.from_root(
    file,
    entry_stop=1,
    ##metadata={"dataset": "DoubleMuon"},
    #schemaclass=BaseSchema,
    schemaclass= NanoAODSchema,
    treepath="/events",
).events()
'''
#events = processor.LazyDataFrame(file['events'], entrystop=20000000)
#out = p.process(events)

fileset = {'JetandStuff': [filename]}
print("fileset")
print(fileset)
'''
output = processor.run_uproot_job(fileset,
    treename='events',
    processor_instance=MyProcessor(),
    executor=processor.futures_executor,
    #executor_args={'workers':opt.cpu},
    executor_args={'schema': None},#processor.LazyDataFrame},
    maxchunks =3,
    chunksize = 1000,
)
'''

output = processor.run_uproot_job(fileset,
    treename='events',
    processor_instance=MyProcessor(),
    executor=processor.futures_executor,
    #executor_args={'workers':opt.cpu},
    executor_args={'schema': None, 'workers':1,},#processor.LazyDataFrame},
    maxchunks =1,
    chunksize = 5000,
)


#print(output)
#print(len(output['nConsti'].value))
#Writeh5(output, 'nConsti_b', 'h5_output')
####Writeh5(output, 'Zuds_weighted_sparse_sshort', 'h5_output')
#print(len(out['ML']['pid']))
#print(len(out['ML']['global']))

print('REMEMBER LOOP MASK')

print(len(output['pid'].value))
print(len(output['data'].value))
print(len(output['coords'].value))
print(output['shape'].value)
print(output['fill_value'].value)
print('COORDS MUST BE TRANSPOSED IN H5 FILE!!!')












