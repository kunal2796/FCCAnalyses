import numpy as np
import h5py
import matplotlib.pyplot as plt

#f = h5py.File('uproot_file.h5', 'r')
#f = h5py.File('h5_Zjets_2104_cut/QG_comb.h5', 'r')
f = np.load('../ConvNet/eval_filesplitter.npz')

pid = f['arr_1']
CNN = f['arr_0'][:,2]
p4 = f['arr_2']
print(np.unique(pid))
#print(pid)
#print(CNN)
#print(p4)

# Here I estimate the luminosity... I should get a better value for the cross section
cross_section = 13000
cross_section_norm = 18616.5
BR = 2*0.156+0.116

def compute_lumi(N):
    return N/cross_section_norm
    #return N/(cross_section*BR)


#quark_mask = (pid==2)|(pid==1)|(pid==3)
#gluon_mask = (pid==21)
'''
strange_mask = (pid==2)
#others_mask = (pid!=2)
down_mask = (pid==0)

p4_s = p4[strange_mask]
#p4_o = p4[others_mask]
p4_d = p4[down_mask]
'''

#1. compute invariant mass -> down select invariant mass via event mask, and down select pid via event mask (pid is double the length so consider only every second element, they are duplicates) -> cut the result into u,d,s event and plot...
# The plan for the Zpeak is to compute the Zpeak for doubly s-tagged events. To this end we should first make a mask for s tagged events.
# For now I require a classifier score of 0.55, though this is just a guess based on the classifier distribution.


print('REMEMBER TO CHANGE THRESHOLDS')

E_bins = np.linspace(0, 50, 25)

jet_mask = (CNN>0.54473346)
pid_cut = pid[jet_mask]

strange_mask = (pid==2)
down_mask = (pid==0)
up_mask = (pid==1)

strange_mask_cut = (pid==2)&(jet_mask)
down_mask_cut = (pid==0)&(jet_mask)
up_mask_cut = (pid==1)&(jet_mask)

E_s = np.nan_to_num(np.histogram(p4[:,0][strange_mask_cut], bins=E_bins)[0]/np.histogram(p4[:,0][strange_mask], bins=E_bins)[0])

Lumi = compute_lumi(len(jet_mask))
print('E_s')
print(E_s)
print(np.histogram(p4[:,0][strange_mask_cut], bins=E_bins)[0])
print(np.histogram(p4[:,0][strange_mask], bins=E_bins)[0])

#bins = np.linspace(90.95,91.50,50)
#bins = np.linspace(0.2,0.75,20)
# Have a look into color scheme at some point
#plt.hist([inv_mass_s, inv_mass_d, inv_mass_u], histtype='step', label=['strange jets','down jets','up jets'], density=False, bins=bins, stacked=True)
plt.plot(E_bins[:-1], E_s, label=r'strange jet efficiency ($\mathbf{\varepsilon_{sig}}$)')
plt.plot(1000, 0, label=r'$\int \mathcal{L}$ = '+str(np.around(Lumi, 3))+r' pb$^{-1}$', linestyle='')
plt.plot(1000, 0, label=r'LodeNet$_s$ $> 0.5447$ single jet', linestyle='')
#plt.hist([inv_mass_s, inv_mass_o], facecolor='b', edgecolor='b', histtype='step', label='strange jets', density=False, bins=bins, stacked=True)
#plt.hist(inv_mass_o, facecolor='r', edgecolor='r', histtype='step', label='down+up jets', density=False, bins=bins, stacked=True)
#plt.yscale('log')
plt.xlim(min(E_bins), max(E_bins))
plt.ylabel(r'Efficiency ($\mathbf{\varepsilon_{sig}})$')
plt.xlabel(r'Jet Energy [GeV]')
#plt.title('Jet {} Distribution'.format('CNN'))
plt.title(r'$\mathbf{FCCee}$ Delphes Sim. - (s-tagged strange jet efficiency $\mathbf{\varepsilon_{sig}}$), $\sqrt{s}$ = 91 GeV')

handles, labels = plt.gca().get_legend_handles_labels()
order = [1,2,0]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc='upper left', fontsize=9)

#plt.legend(loc='upper left', fontsize=9)
plt.savefig('../plots/efficiencyVE.pdf')

