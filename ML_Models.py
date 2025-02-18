#!/usr/bin/env python
# coding: utf-8

# In[60]:


import pandas as pd
import seaborn as sns
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import matplotlib.pyplot as plt


# In[16]:


slow_vars = pd.read_csv("SlowVaryingAmplitudes.csv", header = None)


# In[19]:


slow_vars


# In[20]:


slow_var_arrays = np.asarray(slow_vars)
print(slow_var_arrays.shape)
print(len(slow_var_arrays))


# In[69]:


fig, axes = plt.subplots(3,3, figsize=[15,5], dpi=300, sharex = True)
fig.suptitle("Slow varying amplitudes", fontsize=15)
for i in range(len(slow_var_arrays)):
    j = i//3
    k = i%3
    axes[j,k].plot(slow_var_arrays[i,:], color='r', linewidth=2)
    axes[j,k].set_title(f'a{i}')

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.88, 
                    wspace=0.4, 
                    hspace=0.4)    
plt.savefig("Slow_varying_amplitudes.png")


# In[82]:


#### MLP model ###
'''
This will be used to map the slow varying amplitudes to the surface elevation observations.
Thus, Y(output) in this case will be the CDIP observations.
And, X(inputs) will be the different slow varying amplitudes at a particular time instant.
The number of observations give the number of training examples(N).

We have the matrix of slow varying amplitudes at the different time instants.
'''

X_req = np.transpose(slow_var_arrays)
Y_req = np.asarray(pd.read_csv('z_disp.txt'))[:168]

data = np.append(X_req, Y_req, axis=1)

scaler = MinMaxScaler()
standardized_data = scaler.fit_transform(data)
x_transformed = standardized_data[:,:-1]
y_transformed = standardized_data[:,-1]

x_train_val, x_test, y_train_val, y_test = train_test_split(x_transformed, y_transformed, random_state=0, test_size=0.2)
x_train, x_val, y_train, y_val = train_test_split(x_train_val, y_train_val, random_state=0, test_size=0.3)


# In[83]:


print(x_train.shape)
print(x_val.shape)
print(x_test.shape)


# In[89]:


from tensorflow import keras
from keras.models import Sequential
from keras.layers import Dense


# In[97]:


from tensorflow.keras.optimizers import Adam
from tensorflow.keras.losses import mse


# In[115]:


model = Sequential([tf.keras.Input(shape=(None, x_train.shape[1])),
                    Dense(30, activation = 'relu'),
                    Dense(10, activation = 'relu'),
                    Dense(1)])

model.compile(loss='mean_squared_error',
             optimizer='Adam',
             metrics = ['mean_squared_error'])

history = model.fit(x_train, y_train,
                   epochs = 50,
                   validation_data = (x_val, y_val),
                   verbose=0)


# In[116]:


model.summary()


# In[117]:


train_loss = history.history['loss']
val_loss = history.history['val_loss']
plt.figure(figsize=[6,4])
plt.plot(train_loss, color='r', linewidth=2)
plt.plot(val_loss, color='b', linewidth=2)
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend(['Training loss', 'Validation loss'])


# In[118]:


y_pred = model.predict(x_test).ravel()


# In[121]:


df = pd.DataFrame({'y_true': y_test, 
              'y_preds' : y_pred})

plt.figure(figsize=[8,8])
sns.regplot(data=df, x='y_true',y='y_preds')
plt.savefig('Scatterd_plot_y')

