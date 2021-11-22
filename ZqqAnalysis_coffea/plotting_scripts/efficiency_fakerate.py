import numpy as np
import h5py
import matplotlib.pyplot as plt

#f = h5py.File('uproot_file.h5', 'r')
#f = h5py.File('h5_Zjets_2104_cut/QG_comb.h5', 'r')
f = np.load('../ConvNet/ROC_params.npz')


tpr_s = f['arr_2']
fpr_s = f['arr_5']
thresholds_s = f['arr_8']

def compute_threshold(fpr_s, thresholds_s, fake_rate):
    return thresholds_s[np.argmin(np.abs(fpr_s-fake_rate))]
for i in [0.01, 0.05, 0.1]:
    print('Class Score = '+str(compute_threshold(fpr_s, thresholds_s, i)))
    print('Fake Rate = '+str(fpr_s[thresholds_s==compute_threshold(fpr_s, thresholds_s, i)]))
    print('Efficiency = '+str(tpr_s[thresholds_s==compute_threshold(fpr_s, thresholds_s, i)]))
    print('------------------------------------------------------')

bins = np.linspace(0.05,0.75,20)
#bins = np.linspace(0.2,0.75,20)
print(len(thresholds_s))
#plt.plot(thresholds_s, tpr_s, color='b', label='strange jet efficiency', linestyle="", marker="o")
plt.plot(thresholds_s, tpr_s, color='b', label=r'strange jet efficiency ($\mathbf{\varepsilon_{sig}}$)')
plt.plot(thresholds_s, fpr_s, color='r', label=r'strange jet fake rate ($\mathbf{\varepsilon_{bkg}}$)')
##plt.plot(thresholds_s, tpr_s/fpr_s, color='b', label=r'strange jet efficiency ($\mathbf{\varepsilon_{sig}}$)')
#plt.hist(CNN_o, facecolor='r', edgecolor='r', histtype='step', label='down+up jets', density=False, bins=bins)
#plt.yscale('log')
#plt.xlim([0.086,0.087])
plt.xlim([0,1])
plt.ylabel(r'Efficiency ($\mathbf{\varepsilon}$)')
plt.xlabel(r'Classifier Score ($\mathit{strange}$ $\mathit{node}$)')
plt.title('Jet {} Distribution'.format('CNN'))
plt.title(r'$\mathbf{FCCee}$ Delphes Sim. - (LodeNet $\mathbf{Zuds}$ $\mathbf{\varepsilon_{sig}}$, $\mathbf{\varepsilon_{bkg}}$ ), $\sqrt{s}$ = 91 GeV')
plt.legend(loc='upper right')
plt.savefig('../plots/efficiency_fakerate.pdf')
#plt.savefig('../plots/efficiencyOfakerate.pdf')

