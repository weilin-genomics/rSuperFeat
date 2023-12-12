import pandas as pd
import numpy as np
from scipy import sparse
import pickle
from keras import models,regularizers,optimizers
from keras.layers import Dense
from keras.models import load_model 
from tensorflow import keras
from keras.utils import to_categorical
from sklearn.utils import shuffle,class_weight
from sklearn.preprocessing import binarize
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import collections

# choose gated_genes.hs.csv for human, gated_genes.mm.csv for mouse
genes_used = pd.read_csv('gated_genes.hs.csv', usecols = [0], header = None)
genes_used = list(genes_used.iloc[:, 0])

# load the training data saved by getMatrix function
train_data_file = 'data/train_data.csv'
train_info_file = 'data/train_info.csv'
train_data = pd.read_csv(train_data_file, header = [0], index_col = [0])
train_info = pd.read_csv(train_info_file, header = [0])
train_types = list(train_info.loc[:, 'group'])

x_train = np.array(train_data).T
x_train,train_types,train_names = shuffle(x_train,train_types,train_data.columns,random_state = 0)
label_dict = dict(state0=0, state1=1)
y_train_label = np.array([label_dict[i] for i in train_types])
y_train = to_categorical(y_train_label)

my_class_weight = class_weight.compute_class_weight('balanced', classes=np.unique(y_train_label),y=y_train_label)
my_class_weight = dict(enumerate(my_class_weight))

def plot_acc(history):
    plt.scatter([i for i in range(1, len(history["accuracy"]) + 1)], history["accuracy"], label = "training_accuracy")
    plt.plot([i for i in range(1, len(history["val_accuracy"]) + 1)], history["val_accuracy"], label = "validation_accuracy")
    plt.legend()
    plt.xlabel("Epochs")
    plt.ylabel("Accuracy")
def plot_loss(history):
    plt.scatter([i for i in range(1, len(history["loss"]) + 1)], history["loss"], label = "training_loss")
    plt.plot([i for i in range(1, len(history["val_loss"]) + 1)], history["val_loss"], label = "validation_loss")
    plt.legend()
    plt.xlabel("Epochs")
    plt.ylabel("Loss")
def build_model(input_num_genes,cell_type_num):
    model = models.Sequential()
    model.add(Dense(1,input_shape=(input_num_genes,),activation="relu",kernel_regularizer=regularizers.l1(0.01)))
    model.add(Dense(cell_type_num, activation="softmax")) 
    model.compile(optimizer='Adam',loss='categorical_crossentropy',metrics=["accuracy"])
    return model

# training ...
model=build_model(19202,len(my_class_weight))
history = model.fit(x_train, y_train, batch_size = 128, validation_split = 0.3, shuffle = True, 
                class_weight = my_class_weight, epochs = 30)

# check if stable performance is achieved, or else increase the epoch
plt.figure(figsize = (4,3))
plot_acc(history.history)
plot_loss(history.history)

# check if state1 has positive weights
w2 = pd.DataFrame(model.get_weights()[2])
plt.figure(figsize=(2,0.5))
sns.heatmap(w2)

# check the confusion matrix
y_train_pred = model.predict(x_train, verbose = 0)
y_train_pred = np.argmax(y_train_pred,axis=-1)
cm_train = confusion_matrix(y_train_label, y_train_pred)
np.set_printoptions(precision=2)
cm_train_normalized = cm_train.astype('float')/cm_train.sum(axis=1)[:, np.newaxis]
plt.figure(figsize = (3,1.2))
plt.title('Training Dataset')
sns.heatmap(cm_train_normalized)

# save the model
w1 = pd.DataFrame({'weight':model.get_weights()[0].reshape(-1)}, index = genes_used)
b1 = pd.DataFrame(model.get_weights()[1])
w1.to_csv('models/w1.csv', header=False, index=False)
b1.to_csv('models/b1.csv', header=False, index=False)
model.save('models/model.h5')
