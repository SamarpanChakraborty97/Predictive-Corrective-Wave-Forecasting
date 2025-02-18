import numpy as np
import torch
import os
        
def normal_std(x):
    return x.std() * np.sqrt((len(x) - 1.)/(len(x)))

class Scaler():
    '''
    Performs operations on scaling and transforming the given data
    '''
    def __init__(self, mean, std):
        self.mean = mean
        self.std = std

    def transform(self, data):
        return (data-self.mean)/self.std

    def inverse_transform(self, data):
        return torch.Tensor(data*self.std + self.mean)


class MinMaxScaler_Custom():
    '''
    Performs operations on scaling and transforming the given data
    '''
    def __init__(self, min_val, max_val):
        self.min_val = min_val
        self.max_val = max_val

    def transform(self, data):
        scaled_matrix = np.empty_like(data, dtype=float)
        for i in range(data.shape[1]):  # Iterate over each feature/column
            scaled_matrix[:, i] = 2 * (data[:, i] - self.min_val[i]) / (self.max_val[i] - self.min_val[i]) - 1
        return scaled_matrix

    def inverse_transform(self, data):
        original_matrix = np.empty_like(data, dtype=float)
        for i in range(data.shape[1]):  # Iterate over each feature/column
            original_matrix[:, i] = ((data[:, i] + 1) / 2) * (self.max_val[i] - self.min_val[i]) + self.min_val[i]
        return torch.Tensor(original_matrix)
        

class DataLoaderMulti(object):
    def __init__(self, file_name, train, device, seq_in_len, seq_out_len):
        self.seq_in_len = seq_in_len
        self.seq_out_len = seq_out_len
        self.device = device
        fin = open(file_name)
        self.rawdat = np.loadtxt(fin, delimiter=',').T
        
        self.mean, self.std = self.mean_std(self.rawdat)
        # self.scaler = Scaler(self.mean, self.std)
        # self.dat = self.scaler.transform(self.rawdat)
        # print(self.mean)
        
        self.min, self.max = self.min_max(self.rawdat)
        self.minmaxScaler = MinMaxScaler_Custom(self.min, self.max)
        self.dat = self.minmaxScaler.transform(self.rawdat)
        
        self.n, self.m = self.dat.shape      

        self._split(int(train * self.n))
        
        tmp = self.valid[1] * torch.from_numpy(self.mean).float().expand(self.valid[1].size(0), self.valid[1].size(1), self.m)
        
        self.rse = normal_std(tmp)
        self.rae = torch.mean(torch.abs(tmp - torch.mean(tmp)))

        print(f"The dataset shape is: {self.dat.shape}")
        print(f"The input sequence length is: {self.seq_in_len}")
        print(f"The output sequence length is: {self.seq_out_len}")

        print(f"The train X-dataset shape is: {self.train[0].shape}")
        print(f"The validation X-dataset shape is: {self.valid[0].shape}\n")
        print(f"The train y-dataset shape is: {self.train[1].shape}")
        print(f"The validation y-dataset shape is: {self.valid[1].shape}\n")

    def _split(self, train):

        train_set = range(0, train)
        valid_set = range(train, self.n)
        
        self.train = self._batchify_multi(train_set)
        self.valid = self._batchify_multi(valid_set)

    def mean_std(self, data):
        '''
        Calculates the mean and std of the given data
        '''
        # print(data.shape)
        mean = data.mean(axis=0, keepdims = True)
        # print(mean.shape)
        std = data.std(axis=0, keepdims = True)
        print(std.shape)
        return mean, std

    def min_max(self, data):
        '''
        Calculates the min and max of the given data
        '''
        min_val = np.min(data, axis=0)
        max_val = np.max(data, axis=0)
        
        return min_val, max_val

    def _batchify_multi(self, idx_set):
        n = len(idx_set)-self.seq_in_len-self.seq_out_len+1
        # print(f"Length of the data set: {len(idx_set)}")
        # print(f"In sequence length: {self.seq_in_len}")
        # print(f"Out sequence length: {self.seq_out_len}")
        # print(f"Size of the set: {n}")
        
        X = torch.zeros((n-1, self.seq_in_len, self.m))
        Y = torch.zeros((n-1, self.seq_out_len, self.m))
        for i in range(n-1):
            # print(f"i:{i}")
            start_input = idx_set[i]
            end_input = idx_set[i+self.seq_in_len]
            
            start_output = idx_set[i+self.seq_in_len]
            # print(f"Index of the data in the dataset:{i+self.seq_in_len+self.seq_out_len}\n")
            end_output = idx_set[i+self.seq_in_len+self.seq_out_len]

            # if (i==0 or i==n-2):
            #     print(start_input)
            #     print(end_input)
            #     print(start_output)
            #     print(end_output)
            
            X[i] = torch.from_numpy(self.dat[start_input:end_input, :])
            Y[i] = torch.from_numpy(self.dat[start_output:end_output, :])
        return [X, Y, self.mean, self.std]

    def get_batches(self, inputs, targets, batch_size, shuffle=True):
        length = len(inputs)
        if shuffle:
            index = torch.randperm(length)
        else:
            index = torch.LongTensor(range(length))
        start_idx = 0
        while (start_idx < length):
            end_idx = min(length, start_idx + batch_size)
            excerpt = index[start_idx:end_idx]
            X = inputs[excerpt]
            Y = targets[excerpt]
            
            X = torch.Tensor(X)
            Y = torch.Tensor(Y)
             
            X = X.to(self.device)
            Y = Y.to(self.device)
            
#             yield torch.Tensor(X), torch.Tensor(Y)
            yield X, Y

            start_idx += batch_size
