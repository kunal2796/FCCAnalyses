import matplotlib.pyplot as plt
import numpy as np

d = np.load('../ConvNet/training_params.npz')
train_events = d['arr_4']
loss = d['arr_0']
acc = d['arr_1']
val_loss = d['arr_2']
val_acc = d['arr_3']

# Note that the 6 here is supposed to be the highest epoch to which it trains (which we could get from the training trivially by replacing arr_4)
xaxis = np.linspace(0,6,len(acc))
xaxis_val = np.linspace(0,6,len(val_acc))
plt.figure(figsize=(12, 8))
plt.plot(xaxis, acc, label='Training Accuracy')
plt.plot(xaxis_val, val_acc, label='Validation Accuracy')
plt.ylim([0, 1.1])
#plt.plot([initial_epochs-1,initial_epochs-1], plt.ylim(), label='Start Fine Tuning')
plt.title('Training and Validation Accuracy and Loss')
plt.plot(xaxis, loss, label='Training Loss')
plt.plot(xaxis_val, val_loss, label='Validation Loss')
plt.legend(loc='lower left')
plt.xlabel('epoch')
plt.grid(True)
plt.savefig('../plots/training_curve.pdf')
