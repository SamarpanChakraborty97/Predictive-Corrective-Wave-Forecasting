import numpy as np
import matplotlib.pyplot as plt
import datetime
import matplotlib.dates as mpldts
import calendar
import math
import netCDF4
from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.losses import mse
import tensorflow as tf
import scipy as sp
from sklearn.preprocessing import StandardScaler, MinMaxScaler
# from scikeras.wrappers import KerasRegressor
# from tensorflow.keras.wrappers.scikit_learn import KerasRegressor
import sklearn.metrics
import os
import seaborn as sns
import pandas as pd
from data_prep import *
from scipy.linalg import sqrtm

def buoy_observations(buoy_number_str, time, increment):
    
    stn = buoy_number_str
    dataset = 'realtime' # Enter 'archive' or 'realtime'
    deploy = '34' # If archive dataset, set deployment number from .nc file
    
    if time[:2] == '24':
        time_mod = '00'+':'+time[-2:]
        start_date = '03/24/2024' + ' ' + time_mod # MM/DD/YYYY HH:MM
    else:
        start_date = '03/23/2024' + ' ' + time # MM/DD/YYYY HH:MM 

#     start_date = '03/23/2024 23:35' # MM/DD/YYYY HH:MM
    duration  = increment # Set length of timeseries (minutes)

    qc_level = 2 # Filter data with qc flags above this number 
    
    # Archive
    data_url = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/' + stn + 'p1/' + stn + 'p1_d' + deploy + '.nc'
    # Realtime
    if dataset == 'realtime':
        data_url = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/' + stn + 'p1_xy.nc'
        
    nc = netCDF4.Dataset(data_url)
    # Turn off auto masking
    nc.set_auto_mask(False)
    
    xdisp = nc.variables['xyzXDisplacement'] # Make a numpy array of three directional displacement variables (x, y, z)
    ydisp = nc.variables['xyzYDisplacement']
    zdisp = nc.variables['xyzZDisplacement']
    qc_flag = nc.variables['xyzFlagPrimary']
    filter_delay = nc.variables['xyzFilterDelay']
    start_time = nc.variables['xyzStartTime'][:] # Variable that gives start time for buoy data collection
    sample_rate = nc.variables['xyzSampleRate'][:] # Variable that gives rate (frequency, Hz) of sampling
    end_time = start_time + (len(xdisp)/sample_rate) # Calulate end time for buoy data collection

    # Get station name and number for plot title
    station_name = nc.variables['metaStationName'][:]
    station_title = station_name.tobytes().decode().split('\x00',1)[0]
    
    # Convert to unix timestamp
    def get_unix_timestamp(human_time,dateFormat):
        unix_timestamp = int(calendar.timegm(datetime.datetime.strptime(human_time, dateFormat).timetuple()))
        return unix_timestamp

    # Convert to human readable timestamp
    def get_human_timestamp(unix_timestamp, dateFormat):
        human_timestamp = datetime.datetime.utcfromtimestamp(int(unix_timestamp)).strftime(dateFormat)
        return human_timestamp
    
    data_start = get_human_timestamp(start_time - filter_delay[0],"%m/%d/%Y %H:%M:%S")
    data_end = get_human_timestamp(end_time - filter_delay[0],"%m/%d/%Y %H:%M:%S")
#     print("data_start: " + data_start)
#     print("  data_end: " + data_end)
    
    #Find UNIX timestamps for user human-formatted start/end dates
    unix_start = get_unix_timestamp(start_date,"%m/%d/%Y %H:%M") 
    unix_end = unix_start + (duration * 60) # Create UNIX end stamp by adding duration to 'unix_start'

    # Create specialized array using UNIX Start and End times minus Filter Delay, and Sampling Period (1/sample_rate) 
    # to calculate sub-second time values that correspond to Z-Displacement sampling values
    sample_time = np.arange((start_time - filter_delay[0]), end_time - filter_delay[0],(1/(sample_rate)))

    # Find corresponding start/end date index numbers in 'sample_time' array    
    start_index = sample_time.searchsorted(unix_start) 
    end_index = sample_time.searchsorted(unix_end)
    
    sample_time_cut = sample_time[start_index:end_index]
    
    sample_time_cut *= 1000
    sample_t_cut_ms = sample_time_cut.astype('datetime64[ms]').astype(datetime.datetime)
    
    z = zdisp[start_index:end_index]
    qc = qc_flag[start_index:end_index]

    # Filter out by quality control level
    z = np.ma.masked_where(qc>qc_level,z)
    
    return z
    

def OptimalInterpolation(amp_model_array, y_elev_array, H):
    # print(y_elev_array.shape)
    # print(amp_model_array.shape)
    x_dim = amp_model_array.shape[0]
    # y_dim = y_elev_array.shape[0]
    y_dim = 1
    Pb = createUncertainty(0,1.0,x_dim)
    Ro = createUncertainty(0,1.0,y_dim)
    I = np.eye(x_dim)
    K = computeKalmanGain(Pb,H,Ro)
    amp_corrected_array = amp_model_array + np.matmul(K, (y_elev_array - np.matmul(H,amp_model_array)))
    
    return amp_corrected_array


def createUncertainty(mean, sigma, dim):
    matrix = np.zeros((dim, dim))
    for i in range(dim):
        matrix[i,i] = mean + sigma**2
    
    return matrix


def computeKalmanGain(Pb,H,Ro):
    prod_matrix = np.matmul(np.matmul(H,Pb),H.T) + Ro
    # print(f"Prod matrix is: {prod_matrix}")
    # print(f"H transpose is: {H.T}\n")
    # print(f"Product of Pb and Ht is: {np.matmul(Pb,H.T)}\n")
    K = np.matmul(np.matmul(Pb,H.T),np.linalg.pinv(prod_matrix))
    
    return K

def perturbedObsKF(x_var_file, y_var_file, x_shape_file, y_shape_file):
    x_data = np.loadtxt(x_var_file, delimiter=',')
    x_data_shape = tuple(np.loadtxt(x_shape_file, delimiter=',', dtype=int))

    y_data = np.loadtxt(y_var_file, delimiter=',')
    y_data_shape = tuple(np.loadtxt(y_shape_file, delimiter=',', dtype=int))

    X = x_data.reshape(x_data_shape)
    Y = y_data.reshape(y_data_shape)

    # print(X.shape)
    # print(Y.shape)

    num_time_entries = X.shape[1]
    ensemble_members = X.shape[0]
    num_X_features = X.shape[2]
    num_Y_features = Y.shape[2]

    Ro = createUncertainty(0,1.0,num_Y_features)

    K = np.zeros((num_time_entries, num_X_features, num_Y_features))
    
    for i in range(num_time_entries):
        x_ensemble = X[:,i,:].T
        # print(f"x_ensemble shape:{x_ensemble.shape}")
        x_mean = np.mean(x_ensemble, axis=1).reshape(-1,1)
        x_spread = x_ensemble - x_mean
        x_covariance = (1/(ensemble_members-1)) * np.matmul(x_spread, x_spread.T)

        y_ensemble = Y[:,i,:].T
        y_mean = np.mean(y_ensemble, axis=1).reshape(-1,1)
        y_spread = y_ensemble - y_mean
        y_covariance = (1/(ensemble_members-1)) * np.matmul(x_spread, y_spread.T)

        prod_matrix_2 = (1/(ensemble_members-1)) * np.matmul(y_spread, y_spread.T) + Ro
        # print(prod_matrix_2.shape)
        prod_matrix_1 = (1/(ensemble_members-1)) * np.matmul(x_spread, y_spread.T)
        # print(prod_matrix_1.shape)
        inv_matrix = np.linalg.pinv(prod_matrix_2)
        # print(inv_matrix.shape)
        kalman_gain = np.matmul(prod_matrix_1, inv_matrix)

        K[i] = kalman_gain

    return K


def LETKF(x_var_file, y_var_file, x_shape_file, y_shape_file, obs_file):
    x_data = np.loadtxt(x_var_file, delimiter=',')
    x_data_shape = tuple(np.loadtxt(x_shape_file, delimiter=',', dtype=int))

    y_data = np.loadtxt(y_var_file, delimiter=',')
    y_data_shape = tuple(np.loadtxt(y_shape_file, delimiter=',', dtype=int))

    obs = np.loadtxt(obs_file, delimiter = ',')

    X = x_data.reshape(x_data_shape)
    Y = y_data.reshape(y_data_shape)

    # print(X.shape)
    # print(Y.shape)

    num_time_entries = X.shape[1]
    ensemble_members = X.shape[0]
    num_X_features = X.shape[2]
    num_Y_features = Y.shape[2]

    Ro = createUncertainty(0,1.0,num_Y_features)
    Pb_inv = (ensemble_members - 1) * np.eye(ensemble_members)
    Ro_inv = np.linalg.pinv(Ro)

    Xa_mean = np.zeros((num_time_entries, num_X_features))
    # print(Xa_mean.shape)
    Xa_spread = np.zeros((num_time_entries, num_X_features, ensemble_members))
    
    for i in range(num_time_entries):
        x_ensemble = X[:,i,:].T
        # print(f"x_ensemble shape:{x_ensemble.shape}")
        x_mean = np.mean(x_ensemble, axis=1).reshape(-1,1)
        x_spread = x_ensemble - x_mean
        x_covariance = (1/(ensemble_members-1)) * np.matmul(x_spread, x_spread.T)

        y_ensemble = Y[:,i,:].T
        y_mean = np.mean(y_ensemble, axis=1).reshape(-1,1)
        y_spread = y_ensemble - y_mean
        y_covariance = (1/(ensemble_members-1)) * np.matmul(x_spread, y_spread.T)

        # print(Ro_inv.shape)
        # print(y_spread.shape)
        # print(Pb_inv.shape)

        Pa_inv = Pb_inv + np.matmul(y_spread.T, np.matmul(Ro_inv, y_spread))
        Pa = np.linalg.pinv(Pa_inv)
        Wa = sqrtm((ensemble_members - 1) * Pa)
        prod_matrix_1 = np.matmul(Pa, y_spread.T)
        prod_matrix_2 = np.matmul(Ro_inv, obs[i] - y_mean)
        wa_mean = np.matmul(prod_matrix_1, prod_matrix_2)

        xa_m = x_mean + np.matmul(x_spread, wa_mean)
        # print(xa_m.shape)
        
        Xa_mean[i] = xa_m.reshape(num_X_features)
        Xa_spread[i] = np.matmul(x_spread, Wa)

    return Xa_mean, Xa_spread

def perturbedObservations(obs_file, y_var_file, y_shape_file):
    obs = np.loadtxt(obs_file, delimiter = ',')

    y_data = np.loadtxt(y_var_file, delimiter=',')
    y_data_shape = tuple(np.loadtxt(y_shape_file, delimiter=',', dtype=int))
    Y = y_data.reshape(y_data_shape)

    perturbed_obs = np.zeros((Y.shape[1], Y.shape[0]))

    for i in range(Y.shape[1]):
        perturbed_obs[i,:] = np.random.normal(0,0.1,(1,Y.shape[0]))

    return perturbed_obs
    
            
def reconstructionElevation(amplitudes, omegas, t_present, period):
    
    num_modes = len(omegas)
    # print(f"Nf:{num_modes}")
    
    t_array = np.zeros((amplitudes.shape[0],1))
    elevation_array = np.zeros((amplitudes.shape[0],1))
    sample_rate = 1.28
    dt = 1/sample_rate
    t = t_present
    
    
    for i in range(len(t_array)): 
        rec_elevation = amplitudes[i,0]
        for j in range(num_modes):
            # a = amplitudes[i,j+1] * math.cos(omegas[j] * (t%period))
            a = amplitudes[i,j+1] * math.cos(omegas[j] * (t))
#             print(f"a:{a}")
            # b = amplitudes[i,j+num_modes+1] * math.sin(omegas[j] * (t%period))
            b = amplitudes[i,j+num_modes+1] * math.sin(omegas[j] * (t))
#             print(f"b:{b}")
            rec_elevation += a + b
#         print(t)
        t_array[i] = t
        elevation_array[i] = rec_elevation
        t += 1
    
    t_obs = t_array - t_present
    
    return t_array, t_obs, elevation_array


def create_model(X_req, hidden_layer_sizes=(50, 50), optimizer='adam'):
    model = Sequential()
    model.add(Dense(hidden_layer_sizes[0], input_dim=X_req.shape[1], activation='relu'))
    for units in hidden_layer_sizes[1:]:
        model.add(Dense(units, activation='relu'))
    model.add(Dense(1))
    model.compile(loss='mean_squared_error', optimizer=optimizer)
    return model


def plot_reconstruction_observations(observations, corrections, fontsize, label1, label2, time):
    plt.rcParams['font.size'] = fontsize
    plt.figure(figsize=[15,6], dpi = 225)
    plt.plot(observations,'b-', marker = 'o', mfc = 'b', markersize = 2,label=label1)
    plt.plot(corrections, 'r-', marker = 'o', mfc = 'r', markersize = 2,label=label2)
    plt.grid('minor')
    plt.legend()
    plt.savefig(f"{label1} and {label2} for the latest window at time {time}.png")
    

def plot_complete_duration(observations, corrections, num_mins, increment, fontsize, label1, label2, text):
    plt.rcParams['font.size'] = fontsize
    plt.figure(figsize=[15,6], dpi = 225)
    x_range = min(len(observations), len(corrections))
    X_range_xaxis = np.arange(0, x_range) * (1/(1.28*60))
    plt.plot(X_range_xaxis, observations[:x_range],'b-', marker = 'o', mfc = 'b', markersize = 2,label=label1)
    plt.plot(X_range_xaxis, corrections[:x_range], 'r-', marker = 'o', mfc = 'r', markersize = 2,label=label2)
    plt.xlabel("Minutes forward")
    plt.ylabel("Surface elevation (m)")
    plt.grid('minor')
    plt.xlabel
    plt.legend()
    plt.savefig(f"Comparisons for the entire duration of {num_mins} minutes_increment of {increment} using {text}.png")



def plot_single_values(array, num_mins, increment, fontsize, label):
    plt.rcParams['font.size'] = fontsize
    plt.figure(figsize=[15,6], dpi = 225)
    X_range_xaxis = np.arange(0, len(array)) * (1/(1.28*60))
    plt.plot(X_range_xaxis, array,'b-', marker = 'o', mfc = 'b', markersize = 2,label=label)
    plt.xlabel("Minutes forward")
    plt.ylabel("Errors (m)")
    plt.grid('minor')
    plt.xlabel
    plt.legend()
    plt.savefig(f"Errors in forecast for {num_mins} minutes_increment of {increment}.png")



def time_series_model(train_len, window_len, scaler, amplitudes, time):
    forecast_amps = []
    fig, axes = plt.subplots(11,3, figsize=[35,25], dpi=300, sharex = True)
    fig.suptitle(f"Predictions for a {window_len} steps window", fontsize=40)
    plt.rcParams.update({'font.size': 30})
    amps_scaled = scaler.fit_transform(amplitudes.T)
    
    for i in range(amplitudes.shape[0]):
        file_str = f"best_LSTM_model_amp {i}"
        file_save_name = os.getcwd() + "/" + file_str + "checkpoint.model.keras"
        model_LSTM = keras.models.load_model(file_save_name)
        
        train_timeX = amps_scaled[:,i]
#         x_scaled = scaler.transform(train_timeX.reshape(-1,1))
        x_scaled = train_timeX
        x_model = train_timeX.reshape(1,x_scaled.shape[0], 1)
#         print(x_model.shape)
        x_scaled_list = list(x_model.reshape(x_model.shape[1]))
#         plt.plot(x_scaled)
        
#         y_scaler = MinMaxScaler(feature_range=(-1,1))
#         train_timeY = train_timeX[-output_len:]
#         y_scaled = y_scaler.fit_transform(train_timeY.reshape(-1,1))
#         plt.plot(y_scaled)
   
        forecast_amp = []
        while len(forecast_amp) <= window_len:
            preds_amp = model_LSTM.predict(x_model)
#             print(preds_amp)
#             forecast_amp.extend(y_scaler.inverse_transform(preds_amp[0,:].reshape(-1,1)))
            forecast_amp.extend(preds_amp[0,:])

            x_scaled_list.extend(preds_amp[0,:])
            x_model = np.asarray(x_scaled_list[-train_len:]).reshape(1,train_timeX.shape[0], 1)

        forecast_amps.append(forecast_amp[:window_len])
        
        X_range_xaxis = np.arange(0, len(train_timeX))
        Y_range_xaxis = np.arange(len(train_timeX), len(train_timeX)+window_len)
        j = i//3
        k = i%3
        
        # plt.figure()
        # plt.plot(X_range_xaxis,x_scaled,'r.')
        # plt.plot(Y_range_xaxis,forecast_amp[:window_len],'k.')
        # plt.show()
        
        axes[j,k].plot(X_range_xaxis,train_timeX,'r.')
        axes[j,k].plot(Y_range_xaxis,forecast_amp[:window_len],'k.')
        axes[j,k].set_title(f'a{i}')
    plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.88, 
                        wspace=0.4, 
                        hspace=0.4)    
    plt.savefig(f"Amplitude predictions at {time}_observation_windows.png")
    plt.show()
    forecast_amps_array = scaler.inverse_transform(np.asarray(forecast_amps).T).T
    np.savetxt(f"forecast_amps_next_window_train_len={train_len}",forecast_amps_array,delimiter = ",")
    print("Forecasting done for the next observation window")
    

def time_series_model_whole(train_len, window_len, output_len, scaler, amplitudes, time):
    fig, axes = plt.subplots(11,3, figsize=[35,25], dpi=300, sharex = True)
    fig.suptitle(f"Predictions for a {window_len} steps window", fontsize=40)
    plt.rcParams.update({'font.size': 30})
    
    amps_scaled = scaler.fit_transform(amplitudes.T)

    file_str = "best_model_all_amplitudes_"
    file_save_name = os.getcwd() + "/" + file_str + "checkpoint.model.keras"
    model_LSTM = keras.models.load_model(file_save_name)
        
    train_timeX = amps_scaled 
    x_model = train_timeX.reshape(1,train_timeX.shape[0], 33)
#     x_scaled_list = list(x_model.reshape(x_scaled.shape[1]))
   
    x_scaled_array = train_timeX
    forecast_amps = np.empty((0, train_timeX.shape[1]))
    while len(forecast_amps[:,0]) <= window_len:
        preds_amp = model_LSTM.predict(x_model).reshape(output_len, -1)
        print(preds_amp.shape)
        for j in range(preds_amp.shape[0]):
            forecast_amps = np.vstack((forecast_amps, preds_amp[j,:]))
            x_scaled_array = np.vstack((x_scaled_array, preds_amp[j,:]))
        x_model = np.asarray(x_scaled_array[-train_len:,:]).reshape(1,train_timeX.shape[0], 33)

    forecast_amps = forecast_amps[:window_len,:]
#         forecast_amps.append(forecast_amp[:153])
        
    X_range_xaxis = np.arange(0, len(train_timeX))
    Y_range_xaxis = np.arange(len(train_timeX), len(train_timeX)+window_len)
    for i in range(amplitudes.shape[0]):
        j = i//3
        k = i%3
        axes[j,k].plot(X_range_xaxis,train_timeX[:,i],'r.')
        axes[j,k].plot(Y_range_xaxis,forecast_amps[:,i],'k.')
        axes[j,k].set_title(f'a{i}')   
    plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.88, 
                        wspace=0.4, 
                        hspace=0.4)    
    plt.savefig(f"Amplitude predictions at {time}_observation_windows_whole_lstm.png")
    # plt.show()
    forecast_amps_array = np.asarray(scaler.inverse_transform(forecast_amps).T)
    np.savetxt(f"forecast_amps_next_window_train_len={train_len}",forecast_amps_array,delimiter = ",")
    print("Forecasting done for the next observation window")


def time_series_lstm_whole_latest(device, train_len, window_len, output_len, amplitudes, time_stamp):
    fig, axes = plt.subplots(11,3, figsize=[35,25], dpi=300, sharex = True)
    fig.suptitle(f"Predictions for a {window_len} steps window", fontsize=40)
    plt.rcParams.update({'font.size': 30})

    train_ratio = 0
    data = amplitudes.T
    min_val = np.min(data, axis=0)
    max_val = np.max(data, axis=0)
    minmaxScaler = MinMaxScaler_Custom(min_val, max_val)
    data_scaled = minmaxScaler.transform(data)

    file_str = "best_lstm_model_all_amplitudes_"
    file_save_name = os.getcwd() + "/models" + '/' + file_str + ".pt"
    with open(file_save_name, 'rb') as f:
        model_LSTM = torch.load(f)
    model_LSTM.to(device)
    model_LSTM.eval()

    x_model = torch.Tensor(data_scaled.reshape(1,data_scaled.shape[0], data_scaled.shape[1])).to(device)
   
    x_scaled_array = data_scaled
    forecast_amps = np.empty((0, data_scaled.shape[1]))
    while len(forecast_amps[:,0]) <= window_len:
        preds_amp = model_LSTM(x_model).reshape(output_len, -1).cpu().detach().numpy()
        for j in range(preds_amp.shape[0]):
            forecast_amps = np.vstack((forecast_amps, preds_amp[j,:]))
            x_scaled_array = np.vstack((x_scaled_array, preds_amp[j,:]))
        x_model = torch.Tensor(np.asarray(x_scaled_array[-train_len:,:]).reshape(1,data_scaled.shape[0], 33)).to(device)

    forecast_amps = forecast_amps[:window_len,:]
#         forecast_amps.append(forecast_amp[:153])
        
    X_range_xaxis = np.arange(0, len(data_scaled))
    Y_range_xaxis = np.arange(len(data_scaled), len(data_scaled)+window_len)
    for i in range(amplitudes.shape[0]):
        j = i//3
        k = i%3
        axes[j,k].plot(X_range_xaxis,data_scaled[:,i],'r.')
        axes[j,k].plot(Y_range_xaxis,forecast_amps[:,i],'k.')
        axes[j,k].set_title(f'a{i}')   
    plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.88, 
                        wspace=0.4, 
                        hspace=0.4)    
    plt.savefig(f"Amplitude predictions at {time_stamp}_observation_windows_whole_using_lstm.png")
    plt.close()
    forecast_amps_array = np.asarray(minmaxScaler.inverse_transform(forecast_amps).T)
    np.savetxt(f"forecast_amps_next_window_train_len={train_len}",forecast_amps_array,delimiter = ",")
    print("Forecasting done for the next observation window")


def time_series_gnn_whole_latest(device, train_len, window_len, output_len, amplitudes, time_stamp):
    fig, axes = plt.subplots(11,3, figsize=[35,25], dpi=300, sharex = True)
    fig.suptitle(f"Predictions for a {window_len} steps window", fontsize=40)
    plt.rcParams.update({'font.size': 30})

    train_ratio = 0
    data = amplitudes.T
    min_val = np.min(data, axis=0)
    max_val = np.max(data, axis=0)
    minmaxScaler = MinMaxScaler_Custom(min_val, max_val)
    data_scaled = minmaxScaler.transform(data)

    file_str = "best_gnn_model_all_amplitudes_"
    file_save_name = os.getcwd() + "/models" + '/' + file_str + ".pt"
    with open(file_save_name, 'rb') as f:
        model_GNN = torch.load(f)
    model_GNN.to(device)
    model_GNN.eval()

    x_model = torch.Tensor(data_scaled.reshape(1,data_scaled.shape[0], data_scaled.shape[1])).to(device)
    x_model = torch.unsqueeze(x_model,dim=1)
    x_model = x_model.transpose(2,3)
   
    x_scaled_array = data_scaled
    forecast_amps = np.empty((0, data_scaled.shape[1]))
    while len(forecast_amps[:,0]) <= window_len:
        preds_amp = model_GNN(x_model).reshape(output_len, -1).cpu().detach().numpy()
        for j in range(preds_amp.shape[0]):
            forecast_amps = np.vstack((forecast_amps, preds_amp[j,:]))
            x_scaled_array = np.vstack((x_scaled_array, preds_amp[j,:]))

        
        x_model = torch.Tensor(np.asarray(x_scaled_array[-train_len:,:]).reshape(1,data_scaled.shape[0], 33)).to(device)
        x_model = torch.unsqueeze(x_model,dim=1)
        x_model = x_model.transpose(2,3)

    forecast_amps = forecast_amps[:window_len,:]
#         forecast_amps.append(forecast_amp[:153])
        
    X_range_xaxis = np.arange(0, len(data_scaled))
    Y_range_xaxis = np.arange(len(data_scaled), len(data_scaled)+window_len)
    for i in range(amplitudes.shape[0]):
        j = i//3
        k = i%3
        axes[j,k].plot(X_range_xaxis,data_scaled[:,i],'r.')
        axes[j,k].plot(Y_range_xaxis,forecast_amps[:,i],'k.')
        axes[j,k].set_title(f'a{i}')   
    plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.88, 
                        wspace=0.4, 
                        hspace=0.4)    
    plt.savefig(f"Amplitude predictions at {time_stamp}_observation_windows_whole_using_gnn.png")
    plt.close()
    forecast_amps_array = np.asarray(minmaxScaler.inverse_transform(forecast_amps).T)
    np.savetxt(f"forecast_amps_next_window_train_len={train_len}",forecast_amps_array,delimiter = ",")
    print("Forecasting done for the next observation window")


def time_series_gnn_whole_only(device, train_len, window_len, output_len, amplitudes, time_stamp):
    fig, axes = plt.subplots(11,3, figsize=[35,25], dpi=300, sharex = True)
    fig.suptitle(f"Predictions for a {window_len} steps window", fontsize=40)
    plt.rcParams.update({'font.size': 30})

    train_ratio = 0
    data = amplitudes.T
    min_val = np.min(data, axis=0)
    max_val = np.max(data, axis=0)
    minmaxScaler = MinMaxScaler_Custom(min_val, max_val)
    data_scaled = minmaxScaler.transform(data)

    file_str = "best_gnn_model_all_amplitudes_"
    file_save_name = os.getcwd() + "/models" + '/' + file_str + ".pt"
    with open(file_save_name, 'rb') as f:
        model_GNN = torch.load(f)
    model_GNN.to(device)
    model_GNN.eval()

    x_model = torch.Tensor(data_scaled.reshape(1,data_scaled.shape[0], data_scaled.shape[1])).to(device)
    x_model = torch.unsqueeze(x_model,dim=1)
    x_model = x_model.transpose(2,3)
   
    x_scaled_array = data_scaled
    forecast_amps = np.empty((0, data_scaled.shape[1]))
    while len(forecast_amps[:,0]) <= window_len:
        preds_amp = model_GNN(x_model).reshape(output_len, -1).cpu().detach().numpy()
        for j in range(preds_amp.shape[0]):
            forecast_amps = np.vstack((forecast_amps, preds_amp[j,:]))
            x_scaled_array = np.vstack((x_scaled_array, preds_amp[j,:]))

        
        x_model = torch.Tensor(np.asarray(x_scaled_array[-train_len:,:]).reshape(1,data_scaled.shape[0], 33)).to(device)
        x_model = torch.unsqueeze(x_model,dim=1)
        x_model = x_model.transpose(2,3)

    forecast_amps = forecast_amps[:window_len,:]
#         forecast_amps.append(forecast_amp[:153])
        
    X_range_xaxis = np.arange(0, len(data_scaled))
    Y_range_xaxis = np.arange(len(data_scaled), len(data_scaled)+window_len)
    for i in range(amplitudes.shape[0]):
        j = i//3
        k = i%3
        axes[j,k].plot(X_range_xaxis,data_scaled[:,i],'r.')
        axes[j,k].plot(Y_range_xaxis,forecast_amps[:,i],'k.')
        axes[j,k].set_title(f'a{i}')   
    plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.88, 
                        wspace=0.4, 
                        hspace=0.4)    
    plt.savefig(f"Amplitude predictions at {time_stamp}_observation_windows_whole_using_gnn.png")
    plt.close()
    forecast_amps_array = np.asarray(minmaxScaler.inverse_transform(forecast_amps).T)
    np.savetxt(f"forecast_amps_next_window_train_len={train_len}_forecast_len={window_len}",forecast_amps_array,delimiter = ",")
    print("Forecasting done for the next observation window")
    

def time_series_model_only(train_len, forecast_len, amplitudes, scaler):
    forecast_amps = []
    for i in range(amplitudes.shape[0]):
#         file_str = f"train_len_{train_len}_amplitude_num_{i}_trained_model"
#         print(f"Mode number:{i}")
#         file_str = f"amplitude_num_{i}_trained_model"
        file_str = f"amplitude_model_saves_best_LSTM_train_len_{train_len}_amplitude_num_{i}_trained_model"
#         file_save_name = os.getcwd() + "/amplitude_model_saves2"  + "/best_LSTM_" + file_str + "checkpoint.model.pb"
        file_save_name = os.getcwd() + "/amplitude_model_saves2"  + "/" + file_str + "checkpoint.model.pb"
#         file_save_name = os.getcwd() + "/" + file_str + "checkpoint.model.keras"
        model_LSTM = keras.models.load_model(file_save_name)
        
        train_timeX = amplitudes[i,:]
        x_scaled = scaler.fit_transform(train_timeX.reshape(-1,1))
        
#         plt.figure()
#         plt.plot(train_timeX)
#         plt.show()
        x_model = x_scaled.reshape(1,x_scaled.shape[0], x_scaled.shape[1])
        
        x_scaled_list = list(x_scaled.reshape(x_scaled.shape[0]))
    
        forecast_amp = []
        while len(forecast_amp) <= forecast_len:
            preds_amp = model_LSTM.predict(x_model)
            test_output_preds = scaler.inverse_transform(preds_amp)
            forecast_amp.extend(test_output_preds[0,:])

            x_scaled_list.extend(preds_amp[0,:])
            x_model = np.asarray(x_scaled_list[-train_len:]).reshape(1,x_scaled.shape[0], x_scaled.shape[1])
            
#         X_range_xaxis = np.arange(0, len(train_timeX))
#         Y_range_xaxis = np.arange(len(train_timeX), len(train_timeX)+384)
#         plt.figure()
#         plt.plot(X_range_xaxis,train_timeX,'r')
#         plt.plot(Y_range_xaxis,forecast_amp[:forecast_len],'k')
#         plt.show()
        forecast_amps.append(forecast_amp[:forecast_len])
    forecast_amps_array = np.asarray(forecast_amps)
    np.savetxt(f"forecast_amps_next_window_train_len={train_len}_forecast_len={forecast_len}",forecast_amps_array,delimiter = ",")
    print("Forecasting done for the next observation window")
    
    
def MLP_mapping_results(truth, predictions):
    df = pd.DataFrame({'y_true': truth, 
                       'y_preds' : predictions})
    plt.rcParams["font.size"] = 14
    plt.figure(figsize=[8,6])
    sns.regplot(data=df, x='y_true',y='y_preds')
    ax = plt.gca()
    r,p = sp.stats.pearsonr(df['y_true'], df['y_preds'])
    ax.text(0.05, 0.8, 'r={:.2f}     p={:.2g}'.format(r,p), transform=ax.transAxes)
    print("Predictions mapped to observations for the last few observations in this window.")
    plt.savefig(f'Observations_MLP_mapping_latest_window.png')