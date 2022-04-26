import numpy as np
import sklearn.metrics as metrics
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

f = np.load('../ConvNet/eval_filesplitter.npz')

pid = f['arr_1']
DNN = f['arr_0']
pred_DNN = np.argmax(DNN, axis=1)

cmatrix = metrics.confusion_matrix(pid, pred_DNN, normalize='true') #already the same structure as Vini matrix
cmatrix_col = metrics.confusion_matrix(pid, pred_DNN, normalize='pred')

cmatrix = cmatrix.round(2)
cmatrix_col = cmatrix_col.round(2)

np.savetxt('cmatrix.txt', cmatrix)

print(pid[:20])

DNN_d = DNN[:,0]
pid_d = (np.abs(pid)==0)
DNN_u = DNN[:,1]
pid_u = (np.abs(pid)==1)
DNN_s = DNN[:,2]
pid_s = (np.abs(pid)==2)

print('length of pid_u = '+str(np.sum(pid_u)))
print('length of pid_d = '+str(np.sum(pid_d)))
print('length of pid_s = '+str(np.sum(pid_s)))

fpr_u, tpr_u, thresholds_u = metrics.roc_curve(pid_u, DNN_u)
fpr_d, tpr_d, thresholds_d = metrics.roc_curve(pid_d, DNN_d)
fpr_s, tpr_s, thresholds_s = metrics.roc_curve(pid_s, DNN_s)

#plt.plot(tpr_u, thresholds_u, label='up quarks', linestyle='-', color='b')
plt.plot(fpr_s, thresholds_s, label='fake rate', linestyle='-', color='r')
plt.plot(tpr_s, thresholds_s, label='strange', linestyle='-', color='c')
#plt.plot(1-np.linspace(0, 1, 40), np.linspace(0, 1, 40), linestyle='-', color='k')

#bins = np.linspace(0,1,100)
#plt.hist(tpr_s, facecolor='b', edgecolor='b', histtype='step', label='strange eff', density=False, bins=thresholds_s)
#plt.hist(fpr_s, facecolor='r', edgecolor='r', histtype='step', label='fake rate', density=False, bins=thresholds_s)
plt.xlabel(r'Classifier Score')
plt.ylabel(r'Signal/Background Efficiency ($\mathbf{\varepsilon_{sig}}$/$\mathbf{\varepsilon_{bkg}}$)')
plt.title(r'$\mathbf{FCCee}$ Delphes Sim. - (LodeNet on $\mathbf{Zuds}$ Jets), $\sqrt{s}$ = 91 GeV')
plt.xticks(np.linspace(0, 1, 11))
plt.yticks(np.linspace(0, 1, 11))
plt.grid()
plt.xlim([0,1])
plt.ylim([0,1])
#plt.yscale('log')
plt.legend(loc='best')
plt.savefig('../ConvNet/plots/LodeNet_effandfakes.png')

print("efficiency vs classifier score plot done")