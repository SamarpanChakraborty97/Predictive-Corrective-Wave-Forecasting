import argparse
import math
import time

import torch
import torch.nn as nn
from net import gtnet
import numpy as np
import importlib

from data_prep import *
from trainer import Optim
import os

def evaluate(data, X, Y, model, evaluateL2, evaluateL1, batch_size, mean, std):
    model.eval()
    total_loss = 0
    total_loss_l1 = 0
    n_samples = 0
    predict = None
    test = None

    for X, Y in data.get_batches(X, Y, batch_size, False):
        X = torch.unsqueeze(X,dim=1)
        X = X.transpose(2,3)
        with torch.no_grad():
            output = model(X)
        output = torch.squeeze(output)
        if len(output.shape)==1:
            output = output.unsqueeze(dim=0)
        if predict is None:
            predict = output
            test = Y
        else:
            print(predict.shape)
            print(output.shape)
            predict = torch.cat((predict, output))
            test = torch.cat((test, Y))

        total_loss += evaluateL2(mean + output * std, mean + Y * std).item()
        total_loss_l1 += evaluateL1(mean + output * std, mean + Y * std).item()
        n_samples += (output.size(0) * data.m)

    rse = math.sqrt(total_loss / n_samples) / data.rse
    rae = (total_loss_l1 / n_samples) / data.rae

    predict = (predict).data.cpu().numpy()
    Ytest = test.data.cpu().numpy()
    sigma_p = (predict).std(axis=0)
    sigma_g = (Ytest).std(axis=0)
    mean_p = predict.mean(axis=0)
    mean_g = Ytest.mean(axis=0)
    index = (sigma_g != 0)
    correlation = ((predict - mean_p) * (Ytest - mean_g)).mean(axis=0) / (sigma_p * sigma_g)
    correlation = (correlation[index]).mean()
    return rse, rae, correlation, predict

def train(data, X, Y, model, criterion, lr, wd, batch_size, clip, mean, std, cl, cl_step_size, seq_out_len):
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=wd)
    model.train()
    optimizer.zero_grad()
    total_loss = 0
    n_samples = 0
    task_level = 0
    
    iter = 0
    for X, Y in data.get_batches(X, Y, batch_size, True):
        model.zero_grad()
        X = torch.unsqueeze(X,dim=1)
        X = X.transpose(2,3)
        # Y = torch.unsqueeze(Y,dim=1)
        # Y = Y.transpose(2,3)
        # print(X.shape)
        if iter % params['step_size'] == 0:
            perm = np.random.permutation(range(params['num_nodes']))
        num_sub = int(params['num_nodes'] / params['num_split'])
#         print(f"num_sub:{num_sub}")
#         print(f"perm :{perm}") 

        for j in range(params['num_split']):
            # print(f"j:{j}")
            if j != params['num_split'] - 1:
                id = perm[j * num_sub:(j + 1) * num_sub]
            else:
                id = perm[j * num_sub:]
                
#             print(f"id:{id}")
            id = torch.LongTensor(id).to(device)
            tx = X[:, :, id, :]
            ty = Y[:, :, id]
            output = model(tx,id)
            # print(f"output shape:{output.shape}")
            output = torch.squeeze(output)
            # print(f"output shape:{output.shape}")
            predictions = mean + output * std
            real_values = mean + ty * std

            if iter%cl_step_size==0 and task_level<=seq_out_len:
                task_level +=1
            if cl:
                loss = criterion(predictions[:, :, :task_level], real_values[:, :, :task_level])
            else:
                loss = criterion(predictions, real_values)

            loss.backward()

            torch.nn.utils.clip_grad_norm_(model.parameters(), clip)
            
            total_loss += loss.item()
            n_samples += (output.size(0) * data.m)
            grad_norm = optim.step()

        # if iter%100==0:
        #     print('iter:{:3d} | loss: {:.3f}'.format(iter,loss.item()/(output.size(0) * data.m)))
        iter += 1
        # print(f"The number of batches used are:{iter}")
        
    return total_loss / n_samples

## Parameters for hyperparameter tuning
params= {}

params['data']=os.getcwd() + '/' + 'SlowVaryingAmplitudes_23_30.csv'
params['log_interval']=50
params['optim']='adam'
params['L1Loss']=True
params['device']='cuda'
params['gcn_true']=True
params['buildA_true']=True
params['gcn_depth']=2
params['num_nodes']=33
params['dropout']=0.3
params['subgraph_size']=5
params['node_dim']=40
params['dilation_exponential']=1
params['conv_channels']=16
params['residual_channels']=16
params['skip_channels']=32
params['end_channels']=64
params['in_dim']=1
params['seq_in_len']=100
params['seq_out_len']=50
params['layers']=5
params['batch_size']=16
params['lr']=0.001
params['weight_decay']=0.001
params['clip']=5
params['propalpha']=0.05
params['tanhalpha']=3
params['cl']=False
params['cl_step_size']=50
params['epochs']=250
params['num_split']=1
params['step_size']=10
params['save']=os.getcwd()+'/model'+f"model_layers_{params['layers']}.pt"

device = torch.device(params['device'])
torch.set_num_threads(3)

Data = DataLoaderMulti(params['data'], 0.6, device, params['seq_in_len'], params['seq_out_len'])
mean = Data.train[2]
std = Data.train[3]

model = gtnet(params['gcn_true'], params['buildA_true'], params['gcn_depth'], params['num_nodes'],
                  device, dropout=params['dropout'], subgraph_size=params['subgraph_size'],
                  node_dim=params['node_dim'], dilation_exponential=params['dilation_exponential'],
                  conv_channels=params['conv_channels'], residual_channels=params['residual_channels'],
                  skip_channels=params['skip_channels'], end_channels=params['end_channels'],
                  seq_length=params['seq_in_len'], in_dim=params['in_dim'], out_dim=params['seq_out_len'],
                  layers=params['layers'], propalpha=params['propalpha'], tanhalpha=params['tanhalpha'], layer_norm_affline=False)
model = model.to(device)

print('\nThe recpetive field size is', model.receptive_field)
nParams = sum([p.nelement() for p in model.parameters()])
print('Number of model parameters is', nParams, flush=True)

if params['L1Loss']:
    criterion = nn.L1Loss(reduction="sum").to(device)
else:
    criterion = nn.MSELoss(reduction="sum").to(device)

evaluateL2 = nn.MSELoss(reduction="sum").to(device)
evaluateL1 = nn.L1Loss(reduction="sum").to(device)

best_val = 10000000
optim = Optim(model.parameters(), params['optim'], params['lr'], params['clip'], lr_decay=params['weight_decay'])

print('begin training')
for epoch in range(1, params['epochs'] + 1):
    epoch_start_time = time.time()
    # train_loss = train(Data, Data.train[0], Data.train[1], model, criterion, optim, args.batch_size, args.clip, torch.Tensor(mean).to(device), torch.Tensor(std).to(device))
    train_loss = train(Data, Data.train[0], Data.train[1], model, criterion, params['lr'], params['weight_decay'], params['batch_size'], params['clip'], torch.Tensor(mean).to(device), torch.Tensor(std).to(device), params['cl'], params['cl_step_size'], params['seq_out_len'])
    val_loss, val_rae, val_corr, val_predict = evaluate(Data, Data.valid[0], Data.valid[1], model, evaluateL2, evaluateL1, params['batch_size'], torch.Tensor(mean).to(device), torch.Tensor(std).to(device))
    if epoch % 50 == 0:
        print('| end of epoch {:3d} | time: {:5.2f}s | train_loss {:5.4f} | valid rse {:5.4f} | valid rae {:5.4f} | valid corr  {:5.4f}'.format(epoch, (time.time() - epoch_start_time), train_loss, val_loss, val_rae, val_corr), flush=True)
    # Save the model if the validation loss is the best we've seen so far.

    if val_loss < best_val:
        with open(params['save'], 'wb') as f:
            torch.save(model, f)
        best_val = val_loss


errors = []

vacc = []
vrae = []
vcorr = []

for i in range(10):
    with open(params['save'], 'rb') as f:
        model = torch.load(f)

    vtest_acc, vtest_rae, vtest_corr, val_predictions = evaluate(Data, Data.valid[0], Data.valid[1], model, evaluateL2, evaluateL1, params['batch_size'], torch.Tensor(mean).to(device), torch.Tensor(std).to(device))
    # test_acc, test_rae, test_corr, test_predictions = evaluate(Data, Data.test[0], Data.test[1], model, evaluateL2, evaluateL1, args.batch_size)

    print("validation rse {:5.4f} | validation rae {:5.4f} | validation corr {:5.4f}".format(vtest_acc, vtest_rae, vtest_corr))
    # print("final test rse {:5.4f} | test rae {:5.4f} | test corr {:5.4f}\n".format(test_acc, test_rae, test_corr))
    
    # val_acc, val_rae, val_corr, test_acc, test_rae, test_corr = vtest_acc, vtest_rae, vtest_corr, test_acc, test_rae, test_corr
    val_acc, val_rae, val_corr = vtest_acc, vtest_rae, vtest_corr
    vacc.append(val_acc)
    vrae.append(val_rae)
    vcorr.append(val_corr)

print('\n\n')
print('10 runs average')
print('\n\n')
print("valid\trse\trae\tcorr")
print("mean\t{:5.4f}\t{:5.4f}\t{:5.4f}".format(np.mean(vacc), np.mean(vrae), np.mean(vcorr)))
print("std\t{:5.4f}\t{:5.4f}\t{:5.4f}".format(np.std(vacc), np.std(vrae), np.std(vcorr)))
print('\n\n')

errors.extend([np.mean(vacc), np.mean(vrae), np.mean(vcorr)])
file_save_name = os.getcwd() + "/" + f"num_layers_{params['layers']}.csv"
np.savetxt(file_save_name, errors, delimiter=',')

import matplotlib.pyplot as plt
import random

fig, axes = plt.subplots(11,3, figsize=[35,25], dpi=300, sharex = True)
fig.suptitle("Multi step ahead predictions", fontsize=40)
plt.rcParams["font.family"]="serif"
plt.rcParams['font.size']=30

random_dataset = random.randint(0,Data.valid[0].shape[0])
for i in range(Data.valid[0].shape[2]):
    j = i//3
    k = i%3
    x_range = np.arange(0,len(Data.valid[0][random_dataset,:,i]))
    y_range = np.arange(len(Data.valid[0][random_dataset,:,i]),len(Data.valid[0][random_dataset,:,i])+len(Data.valid[1][random_dataset,:,i]))
    axes[j,k].plot(x_range, Data.valid[0][random_dataset,:,i], 'k-.',linewidth=2, label="x-data")
    axes[j,k].plot(y_range, Data.valid[1][random_dataset,:,i], color='b', marker = 'h', markersize = 3, label="truth")
    axes[j,k].plot(y_range, val_predictions[random_dataset,:,i], color='r', marker = 'h', markersize = 3, label='prediction')
    axes[j,k].set_title(f'a{i}')
    axes[j,k].tick_params(axis='both', which='major', labelsize=18)

    axes[j,k].legend(ncol = 3, fontsize = 14)

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.88, 
                    wspace=0.4, 
                    hspace=0.4)    
plt.savefig(f"multi step ahead predictions_during training_num_layers_{params['layers']}.png")
plt.show()