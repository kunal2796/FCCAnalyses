import numpy as np
import h5py
import matplotlib.pyplot as plt

#f = h5py.File('uproot_file.h5', 'r')
#f = h5py.File('h5_Zjets_2104_cut/QG_comb.h5', 'r')
f = np.load('../ConvNet/ROC_params.npz')


tpr_s = f['arr_2']
fpr_s = f['arr_5']
thresholds_s = f['arr_8']

# Note that I plotted this and it looked a bit weird...
def Asimov_estimate(fpr_s, tpr_s):
    return np.sqrt(2*((fpr_s+tpr_s)*np.log(1+tpr_s/fpr_s)-tpr_s))

def compute_threshold(fpr_s, thresholds_s, fake_rate):
    return thresholds_s[np.argmin(np.abs(fpr_s-fake_rate))]
for i in [0.01, 0.05, 0.1, 0.665]:
    print('Class Score = '+str(compute_threshold(fpr_s, thresholds_s, i)))
    print('Fake Rate = '+str(fpr_s[thresholds_s==compute_threshold(fpr_s, thresholds_s, i)]))
    print('Efficiency = '+str(tpr_s[thresholds_s==compute_threshold(fpr_s, thresholds_s, i)]))
    print('------------------------------------------------------')

bins = np.linspace(0.05,0.75,20)
#bins = np.linspace(0.2,0.75,20)
print(len(thresholds_s))
#plt.plot(thresholds_s, tpr_s, color='b', label='strange jet efficiency', linestyle="", marker="o")
##plt.plot(thresholds_s, tpr_s, color='b', label=r'strange jet efficiency ($\mathbf{\varepsilon_{sig}}$)')
##plt.plot(thresholds_s, fpr_s, color='r', label=r'strange jet fake rate ($\mathbf{\varepsilon_{bkg}}$)')
plt.plot(thresholds_s, tpr_s/np.sqrt(fpr_s), color='b', label=r'sig/bkg efficiency ($\frac{\mathbf{\varepsilon_{sig}}}{\sqrt{\mathbf{\varepsilon_{bkg}}}}$)')
plt.plot(np.repeat(0.665, len(np.linspace(0,11))), np.linspace(0,11), color='k', label='Adopted threshold = 0.665')
#plt.plot(np.repeat(0.67, len(np.linspace(0,11))), np.linspace(0,11), color='k')
#plt.plot(np.repeat(0.68, len(np.linspace(0,11))), np.linspace(0,11), color='k')
#plt.plot(thresholds_s, Asimov_estimate(fpr_s, tpr_s), color='r', label=r'Asimov estimate')
#plt.hist(CNN_o, facecolor='r', edgecolor='r', histtype='step', label='down+up jets', density=False, bins=bins)
#plt.yscale('log')
#plt.xlim([0.086,0.087])
plt.xlim([0,1])
#plt.ylabel('Number of Jets')
plt.xlabel(r'Classifier Score ($\mathit{strange}$ $\mathit{node}$)')
plt.title('Jet {} Distribution'.format('CNN'))
plt.title(r'$\mathbf{FCCee}$ Delphes Sim. - (LodeNet $\mathbf{Zuds}$ $\frac{\mathbf{\varepsilon_{sig}}}{\mathbf{\varepsilon_{bkg}}}$ ), $\sqrt{s}$ = 91 GeV')
plt.legend(loc='upper left')
#plt.savefig('../plots/efficiency_fakerate.pdf')
plt.savefig('../plots/efficiencyOfakerate.pdf')

