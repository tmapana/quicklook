import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os.path
import scipy
import math
import glob
import os

from scipy import signal
from scipy.signal import firwin, lfilter, iirnotch
from scipy.interpolate import interp1d

# CONSTANTS
c = 299792458
BOLTZMANN = 1.380649e-23

# RED PITAYA
RP_CLK = 125e6
RP_ADC_VPP = 2      # 2V peak-to-peak on this model red pitaya
RP_ADC_BITS = 14    # with 14-bit ADC
ADC_BITS_VOLTS = 2*RP_ADC_VPP/math.pow(2, RP_ADC_BITS)
DR_MAX = 10*pow(2, 20)  # 10 MBps
BYTES_PER_WRITE = 4  # IQ 16 bit
N_CHANNELS = 2
N_COUNTER = 75
FD = pow(2, 24) - 1
RF_DIV = 4
DEPRESSION_ANGLE = 20 # angle of depression of antennas
RAMP_BANDWIDTH = 175e6

#-----------------------------------------------------------------------------------------------------#
def next_pow_two(number):
    '''
    Returns the next power of two.
    '''
    return int(pow(2, np.ceil(np.log2(number))))

#-----------------------------------------------------------------------------------------------------#
def linear2db(linear):
    '''
    Conversion from linear to decibels.
    '''
    return 10*np.log10(linear)

#-----------------------------------------------------------------------------------------------------#
def get_timestamp(series, div=1e6):
    return pd.DataFrame({'timestamp': (series.astype(np.int64)//div).astype(np.int64)})

#-----------------------------------------------------------------------------------------------------#
def load_motion_data(root_directory, parameter):
    # open file stored at existing directory and read as binary
    motion_file_name = glob.glob(os.path.join(root_directory, '*.json'))[0]  

    # read motion dataframe from file
    motion_df = pd.read_json(
        motion_file_name, orient='records', lines=True, precise_float=True)
    motion_df = motion_df.assign(
        timestamp=get_timestamp(motion_df['utc'], div=1))
    mean_height = np.mean(
        motion_df['u'].values)  # average height above base unit
    mean_velocity = np.mean(
        motion_df['abs_vel'].values, where=[True])  # average platform velocity
    synthetic_aperture = np.hypot(
        motion_df['e'].values[-1], motion_df['n'].values[-1])

    if parameter == 'height':
        return mean_height
    elif parameter == 'velocity':
        return mean_velocity
    elif parameter == 'aperture':
        return synthetic_aperture
    else:
        print('Unknown parameter option, please confirm.')
        exit()

#-----------------------------------------------------------------------------------------------------#
def get_rf_freq(FN):
    '''
    Returns the RF frequency associated with a given fractional numerator
    '''
    return RP_CLK*(N_COUNTER + FN/FD)/RF_DIV

#-----------------------------------------------------------------------------------------------------#
def fast_time_fft(data, ns_fft):
  print("\nPerforming range_compression...")
  data = np.fft.fftshift(np.fft.fft(data, ns_fft, axis=0), axes=0)

  return data

#-----------------------------------------------------------------------------------------------------#
def window(data, switch_mode, ns_pri):
    '''
    Attempt at windowing the data to remove high sidelobes
    Using Hamming window function
    '''
    window_function = np.hamming(ns_pri)
    # window_function = np.blackman(ns_pri)

    window_function = np.array([window_function])

    if switch_mode  == 1:
      channel_1a = data[:,:,0]
      channel_1b = data[:,:,1]

      channel_1a = np.multiply(channel_1a, np.transpose(window_function))
      channel_1b = np.multiply(channel_1b, np.transpose(window_function))

      data_out = np.dstack((channel_1a, channel_1b))
    
    if switch_mode  == 2:
      channel_2a = data[:,:,0]
      channel_2b = data[:,:,1]

      channel_2a = np.multiply(channel_2a, np.transpose(window_function))
      channel_2b = np.multiply(channel_2b, np.transpose(window_function))

      data_out = np.dstack((channel_2a, channel_2b))
          
    elif switch_mode == 3:
      channel_1a = data[:,:,0]
      channel_1b = data[:,:,1]
      channel_2a = data[:,:,2]
      channel_2b = data[:,:,3]

      channel_1a = np.multiply(channel_1a, np.transpose(window_function))
      channel_1b = np.multiply(channel_1b, np.transpose(window_function))
      channel_2a = np.multiply(channel_2a, np.transpose(window_function))
      channel_2b = np.multiply(channel_2b, np.transpose(window_function))

      data_out = np.dstack((channel_1a,    # 1a
                            channel_1b,    # 1b
                            channel_2a,    # 2a
                            channel_2b))   # 2b

    return data_out

#-----------------------------------------------------------------------------------------------------#
def haversine(lon_A, lat_A, lon_B, lat_B):
  '''
  This function is going to make use of the Haversine formula to calculate the distance between two points
  Process is as follows:
    - convert latitude and longitude to radians
    - calculate radial distance between lon_A and lon_B
    - calculate radial distance between lat_A and lat_B
    - the result of the Haversine formula is calculated by making use of the law of cosine
  '''
  lon1 = math.radians(lon_A)
  lon2 = math.radians(lon_B)
  lat1 = math.radians(lat_A)
  lat2 = math.radians(lat_B)

  dlon = lon2 - lon1
  dlat = lat2 - lat1
  
  a = math.pow(math.sin(dlat/2), 2) + \
    math.cos(lat1)*math.cos(lat2)*math.pow(math.sin(dlon/2), 2)
  
  c = math.atan2(math.sqrt(a), math.sqrt(1-a))

  R = 6371  # Earth's approximate radius in kilometres 

  # return distance by which to compensate in meters
  distance = 0  # default case
  if dlat > 0:    # provide an error margin
    distance = -1*R*c*1000
  if dlat < 0:
    distance = 1*R*c*1000

  # distance = 0  # debug statement
  return distance # in metres

#-----------------------------------------------------------------------------------------------------#
def trimmer(data, min_chunk, max_chunk, min_range_bin, max_range_bin):
  data = data[min_range_bin:max_range_bin, min_chunk:max_chunk]
  
  print('Dataset trimmed.')
  return data

#-----------------------------------------------------------------------------------------------------#
def cfar_filter(data, n_chunks, n_training_cells):
  # step 1: set parameters to use for the CFAR filter
  kernel_dims = [1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1]
  pfa = 1e-6  # probability of false alarm

  N = np.sum(np.sum(kernel_dims))  # number of training cells
  alpha = N*(pow(pfa, -1/N) - 1)  # threshold gain

  signal_power = pow(abs(data), 2) # signal amplitude

  # cfar
  dims = data.shape
  kernel = np.zeros((100,100,))
  kernel[1:len(kernel_dims),:,:] = kernel_dims
  mask = np.fft.fft2(kernel)
  noise = np.multiply(np.conjugate(mask), np.fft.fft2(signal_power))
  noise = np.fft.ifft(noise)

  threshold = noise*alpha
   
  '''threshold = 0
  for rng_line in range(0, data.shape[1]):
    for az_pt in range(0, data.shape[0]):
      temp_pow_lvl = pow(abs(data[az_pt, rng_line]), 2)  # = I^2 + Q^2
      threshold = threshold + temp_pow_lvl
    
  threshold = threshold * (1/n_training_cells)'''
  print("Threshold level:", threshold)

  # pass through the data again and filter given the threshold

  rfi_count = 0
  data_cfar = data
  for rng_line in range(0, data.shape[1]):
    for az_pt in range(0, data.shape[0]):
       pow_lvl = pow(abs(data_cfar[az_pt, rng_line]), 2)

       if (pow_lvl[0] > threshold[0]) and (pow_lvl[1] > threshold[1]):
          data_cfar[az_pt, rng_line] = np.zeros((1,1)) # I think...
          rfi_count += 1

  print("RFI count:", rfi_count)
  return data_cfar
   
#-----------------------------------------------------------------------------------------------------#
def cfar(root_directory, data, ns_fft, switch_mode, n_chunks, ns_cpi, ns_guard=5, ns_train=25, n_sigma=10):
  # list of indices for moving window
  template_indx = np.arange(ns_guard+1, ns_guard+ns_train+1, 1, dtype='int')
  template_indx = np.concatenate([np.flipud(-template_indx), template_indx])

  # number of channels
  n_datasets = 2
  if switch_mode==3:
    n_datasets = 4

  avg_map = np.zeros((ns_fft, ns_cpi, n_datasets))
  threshold = np.zeros((1, ns_cpi, n_datasets))

  # window with boxcar to ensure minimum number of zeros in mask
  # this is a tradeoff between sidelobe level of mask's impulse response
  # and the efficacy of the mask at zeroing interference.
  # unlike usual, not windowing will result in lower sidelobes, since there are fewer zeros in the mask.
  window = np.transpose(np.repeat(np.ones(ns_cpi)[:, np.newaxis], ns_fft, axis=1))
  window = np.repeat(window[:,:,np.newaxis], n_datasets, axis=2)

  s_col = 0  # start col
  e_col = 1  # end col
  avg_map += abs(np.fft.fft(np.multiply(data[:, s_col: e_col, :], window), n=ns_cpi, axis=1))
  del s_col, e_col, window

  # average along range axis
  avg_map = np.mean(avg_map, axis=0, keepdims=True)

  # determine threshold along Doppler axis
  for i in range(0, ns_cpi):
    indx = template_indx + i
    temp_data = np.take(avg_map, indx, axis=1, mode='wrap')
    threshold[0, i, :] = np.mean(temp_data, axis=1) + n_sigma*np.std(temp_data, axis=1)
  del template_indx, temp_data, indx

  fig, ax = plt.subplots(2, 1)
  ax[0].plot(2*linear2db(np.fft.fftshift(avg_map[0, :, 0])), label='Averaged Doppler')
  ax[0].plot(2*linear2db(np.fft.fftshift(threshold[0, :, 0])), label='Threshold')

  # apply threshold and quantise values, forming Doppler mask
  avg_map -= threshold
  avg_map[avg_map > 0] = 0
  avg_map[avg_map < 0] = 1
  del threshold

  # interpolate Doppler mask from CPI to N_CHUNKS
  x = np.linspace(0, ns_cpi, ns_cpi, endpoint=False).astype('int')
  func = interp1d(x, avg_map, kind='linear', axis=1)
  x = np.linspace(0, ns_cpi, n_chunks, endpoint=False).astype('int')
  avg_map = func(x)
  del func

  ax[0].plot(x, np.fft.fftshift(avg_map[0, :, 0]), label='Doppler Mask')

  ax[0].set_xlabel('Doppler Frequency [Hz]')
  ax[0].grid()
  ax[0].legend()

  impulse_response = np.fft.fftshift(2*linear2db(abs(np.fft.ifft(avg_map[0, :, 0]))))
  ax[1].plot(impulse_response - max(impulse_response))
  ax[1].grid()
  plt.tight_layout()
  plt.savefig(os.path.join(root_directory, 'quicklook/averaged_doppler.png'))

  # print percentage of bins blanked
  for i in range(n_datasets):
    loss = 100 - (sum(avg_map[0, :, i])/n_chunks * 100)
    print('Channel', i, ': ' + str(np.round(loss, 2)) + '% blanked')

  # repeat Doppler mask for all range
  avg_map = np.repeat(avg_map, ns_fft, axis=0)

  # apply range-Doppler mask to data
  spectrum = np.fft.fft(data, n=n_chunks, axis=1)
  spectrum *= avg_map
  data = np.fft.ifft(spectrum, n=n_chunks, axis=1)

  return data

#-----------------------------------------------------------------------------------------------------#
def cfar_detector(root_directory, data, ns_fft, n_chunks, n_guard, n_train, pfa):
  '''
  Cell-averaging Constant False Alarm Rate (CFAR) algorithm to detect peaks
  Python implementation of the Matlab default CFAR algorithm
  https://www.mathworks.com/help/phased/ug/constant-false-alarm-rate-cfar-detection.html
  n_guard: number of guard cells
  n_train: number of training cells
  pfa: probability of false alarm
  '''

  # n_cells = n_chunks
  n_cells = ns_fft
  n_guard_half = round(n_guard / 2)
  n_train_half = round(n_train / 2)
  n_side = n_guard_half + n_train_half

  print("\nExecuting CFAR detector...", "\nNum cells:", n_cells, "\nNum train:", n_train, "\nNum guard:", n_guard, "\nPFA:", pfa)

  # threshold factor
  alpha = n_train * (pfa**(-1/n_train) - 1)

  data_peak_cells = []
  data_new = data

  for cell in range(n_side, n_cells-n_side):
    for y in range(0, n_chunks):
      if (cell != (cell-n_side+np.argmax(data[cell-n_side:cell+n_side+1, y, :]))).all():
        # print(cell, np.argmax(data[cell-n_side:cell+n_side+1, y], axis=0)) # test statement
        continue
      # print(cell, np.argmax(data[cell-n_side:cell+n_side+1, y, :])) # test statement
      
      sum1 = (np.sum(data[cell-n_side:cell+n_side+1, y]))
      sum2 = (np.sum(data[cell-n_guard_half:cell+n_guard_half+1, y]))

      # estimate of noise power in training cells
      p_noise = (sum1 - sum2) / n_train

      # estimate of threshold
      threshold = alpha * p_noise

      # CFAR detection
      if (abs(data[cell, y]) > threshold).all():
        data_peak_cells.append(cell)
        data_new[cell, y] = 0
  
  data_peak_cells = np.array(data_peak_cells, dtype=int)
  print("Number of deatected RFI peaks:", data_peak_cells.shape[0])  # test statement

  return data_new

#-----------------------------------------------------------------------------------------------------#
def notch_filter(data, n_chunks, intergerence_frequency, quality_factor, sampling_rate):
  # Notch filter polynomials
  nyquist_rate = 0.5 * sampling_rate
  w0 = intergerence_frequency / nyquist_rate
  b, a = iirnotch(w0=w0, Q=quality_factor, fs=sampling_rate)

  # Frequency response, [frequency, magnitude]
  freq, h = signal.freqz(b, a, fs=sampling_rate)

  # data_filtered = signal.filtfilt(b, a, data)
  h = h.reshape(n_chunks)
  data_filtered = np.dot(h, data)

  return data_filtered

'''def notch_filter(data, sampling_rate, intereference_freq, notch_width):
  nyquist_freq = sampling_rate / 0.5
  normalized_interference_freq = intereference_freq / nyquist_freq
  normalized_notch_width = notch_width / nyquist_freq

  filter_order = 1
  filter_coeffs = firwin(filter_order, [normalized_interference_freq - 0.5*normalized_notch_width,
                                        normalized_interference_freq + 0.5*normalized_notch_width],
                                        window='hamming')

  #  Apply notch filter
  filtered_data = lfilter(filter_coeffs, 1.0, data)

  return filtered_data'''