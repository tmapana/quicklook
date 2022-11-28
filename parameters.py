import pandas as pd
import numpy as np
import os.path
import math
import glob
import os

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
  data = data[min_range_bin:max_range_bin, min_chunk, max_chunk]
  
  
  print('Dataset trimmed.')
  return data
