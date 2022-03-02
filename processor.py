"""
This is a quicklook processor which is used to give a quick
snapshot of miloSAR data for review

Author:   Tlotliso Mapana
Date:     28 September 2021
"""

from cProfile import label
from cmath import phase
from logging import root
import os
from sqlite3 import Timestamp
import sys
import math
import glob
import json
from matplotlib.colors import LogNorm
import numpy as np
import configparser
from numpy.fft import fftfreq
from numpy.lib.function_base import add_newdoc_ufunc
import matplotlib.pyplot as plt
from numpy.linalg import matrix_rank
import pandas as pd
from scipy import interpolate

#TODO: Motion compensation
#TODO: Basic SAR processor

def main():
  '''
  This is the entry point of the processor
  The main function calls these functions:
    1. load_configuration_parameters
    2. load_data
    3. range_doppler
    4. rti
  '''

  # guide user for correct input
  if (len(sys.argv) < 2):
    print("Usage: processor.py dataset_root_directory [OPTIONS]")
    exit(-1)  # exit on incorrect input

  # the first command-line argument is the root directory of the miloSAR data
  root_directory = sys.argv[1]
  
  # open file and read in the data into an array
  if os.path.exists(root_directory):
    # create output directory for processed data
    if not os.path.exists(os.path.join(root_directory, 'quicklook')):
      os.mkdir(os.path.join(root_directory, 'quicklook'))

    # Extract configuration parameters from setup.ini and summary.ini files
    time_stamp, switch_mode, n_seconds, presummed_prf, n_pri, max_range \
       = load_configuration_parameters(root_directory=root_directory)

    x_axis = np.linspace(0, n_seconds, n_pri, endpoint=False)  # time axis
    y_axis = np.linspace(0, max_range, (int)(n_seconds*presummed_prf), endpoint=False)  # range axis

    # load data from file and pack into a data matrix
    data = load_data(root_directory=root_directory, time_stamp=time_stamp, switch_mode=switch_mode, \
                         n_seconds=n_seconds, presummed_prf=presummed_prf, n_pri=n_pri)

    # perform motion compensation
    #data = motion_compensation(data=raw_data, x_axis=x_axis, switch_mode=switch_mode, root_directory=root_directory)

    # OPTIONS
    if len(sys.argv) > 2:
      for opt in range(2, len(sys.argv)):
        if sys.argv[opt] == 'raw':
          # save .raw file
          save_raw(data, root_directory, time_stamp, switch_mode)
        
        elif sys.argv[opt] == 'rd':
          # produce a range-Doppler plot and save to directory
          range_doppler(data, root_directory, time_stamp, switch_mode, \
                        n_seconds, n_pri, presummed_prf, max_range)

        elif sys.argv[opt] == 'rti':
          # produce range-time plot and save to diectory
          rti(data, root_directory, time_stamp, switch_mode, \
              n_seconds, n_pri, presummed_prf, max_range)
    
        elif sys.argv[opt] == 'hist':
          # produce a histogram of the raw data
          plot_histogram(data, root_directory=root_directory)

        else:
          print("\nUnidentified function selected.")
          print("Usage: processor.py dataset_root_directory [OPTIONS]")
          print("  OPTIONS:")
          print("    raw:    Save .raw file for observing data more accurately")
          print("    rd:     Plot range-Doppler map for each polarization")
          print("    rti:    Plot and save range time intensity plot for each polarization")
          print("    hist:   Plot histogram of data\n")

    else:
      print("\nUsage: processor.py dataset_root_directory [OPTIONS]")
      print("  OPTIONS:")
      print("    raw:    Save .raw file for observing data more accurately")
      print("    rd:     Plot range-Doppler map for each polarization")
      print("    rti:    Plot and save range time intensity plot for each polarization")
      print("    hist:   Plot histogram of data\n")
              
  else:
    print("Directory", root_directory, "does not exist!")
    exit(-1)


def load_configuration_parameters(root_directory):
  '''
  Important information about the data is stored in the configuration file of each data directory
  This function reads in the config file and extracts the following experiment parameters:
    1. switch_mode = integer indicating which channel(s) was used during operation
    2. n_seconds = length of data take in seconds
    3. prf = pulse repetition frequency
    4. sampling_rate = sampling rate of ADC
    5. n_pris = factor applied to the prf when writing data to physical disk
    6. start_index = position to start recording data
    7. end_index = position to end recording
  '''
  summary = configparser.ConfigParser()
  summary.read(os.path.join(root_directory, "summary.ini"))

  time_stamp = summary['general']['time_stamp']
  switch_mode = int(summary['dataset']['switch_mode'])
  n_seconds = int(summary['dataset']['n_seconds'])
  prf = int(summary['dataset']['prf'])
  sampling_rate = float(summary['dataset']['sampling_rate'])
  presumming_factor = int(summary['integration']['n_pris'])   # TODO: define clearly what this is
  start_index = int(summary['integration']['start_index'])
  end_index = int(summary['integration']['end_index'])
  
  setup = configparser.ConfigParser()
  setup.read(os.path.join(root_directory, "setup.ini"))
  max_range_string = setup['geometry']['max_range']  # 900.0; [m]
  max_range = int(max_range_string.split('.')[0])

  # realisable prf of data file stored in memory as data rate is halved when writing to memory
  presummed_prf = prf//presumming_factor # TODO rename to more understandable term

  # number of range bins recorded for each pri
  n_pri = end_index - start_index + 1
  
  print("Configuration parameters loaded.")

  return time_stamp, switch_mode, n_seconds, presummed_prf, n_pri, max_range


def load_data(root_directory, time_stamp, switch_mode, n_seconds, presummed_prf, n_pri):
  '''
  This function reads in raw miloSAR data stored in a .bin file in the given directory
  '''
  # open file stored at existing directory and read as binary
  file_name = os.listdir(root_directory)[0]
  with open(os.path.join(root_directory,file_name), 'rb') as f:
    # read data into an array as signed 16-bit integer values
    data_array = np.fromfile(f, np.int16)
  
  '''
  data_array is an array of interleaved IQ data samples.
  Each sample contains I and Q data point from two channels, A and B
  Need to extract each channel data respectively and reshape data into a 2D matrix
  '''
  # each data sample contains [Q_a, I_a, Q_b, I_b]
  channel_a_data = data_array[1::4] + 1j*data_array[0::4]
  channel_b_data = data_array[3::4] + 1j*data_array[2::4]

  # calculate the number of samples recorded per channel
  if switch_mode==3:
    presummed_prf = int(presummed_prf//2)
  n_samples_per_channel = n_seconds*presummed_prf*n_pri

  '''
  The actual file stored in memory has more bytes than what is calculated
  The file size in memory should be:
  n_seconds * presummed_prf * sample_range * 4(number of data streams)
  That means there are extra samples which need to be discarded.
  '''
  # discard extra samples at the end
  channel_a_data = channel_a_data[:n_samples_per_channel]
  channel_b_data = channel_b_data[:n_samples_per_channel]

  # reshape array into a 2D matrix for each channel in order to perform fft
  channel_a_data = np.transpose(channel_a_data.reshape(
    n_seconds*presummed_prf, n_pri))
  channel_b_data = np.transpose(channel_b_data.reshape(
    n_seconds*presummed_prf, n_pri))
  
  # Pack data into a matrix
  data = np.array([])
  if (switch_mode == 1) or (switch_mode == 2):
    # Only two polarizations recorded for either switch mode 1 or 2
    data = np.dstack((channel_a_data, channel_b_data))

  elif (switch_mode == 3):
    # four polarizations recorded
    data = np.dstack((channel_a_data[:, 0::2],  # 1a
                    channel_b_data[:, 0::2],    # 1b
                    channel_a_data[:, 1::2],    # 2a
                    channel_b_data[:, 1::2]))   # 2b

  print("Dataset", time_stamp, "loaded successfully.")
  return data


def save_raw(data, root_directory, time_stamp, switch_mode):
  '''
  Not currently saving raw data as intended
  TODO: Rewrite this function to produce correct results
  '''
  if (switch_mode == 1):
    # only two polarizations recorded, other channels terminated
    channel_1a = data[:,:,0]
    channel_1b = data[:,:,1]
    with open(os.path.join(root_directory, time_stamp + '_1a_i_.raw'), 'wb') as raw_1a_i:
      (channel_1a.real.astype(np.int16)).tofile(raw_1a_i)
    with open(os.path.join(root_directory, time_stamp + '_1a_q_.raw'), 'wb') as raw_1a_q:
      (channel_1a.imag.astype(np.int16)).tofile(raw_1a_q)
    with open(os.path.join(root_directory, time_stamp + '_1b_i_.raw'), 'wb') as raw_1b_i:
      (channel_1b.real.astype(np.int16)).tofile(raw_1b_i)
    with open(os.path.join(root_directory, time_stamp + '_1b_q_.raw'), 'wb') as raw_1b_q:
      (channel_1b.imag.astype(np.int16)).tofile(raw_1b_q)
  
  elif (switch_mode == 2):
    # only two polarizations recorded, other channel terminated
    channel_2a = data[:,:,0]
    channel_2b = data[:,:,1]
    with open(os.path.join(root_directory, time_stamp + '_2a_i_.raw'), 'wb') as raw_2a_i:
      (channel_2a.real.astype(np.int16)).tofile(raw_2a_i)
    with open(os.path.join(root_directory, time_stamp + '_2a_q_.raw'), 'wb') as raw_2a_q:
      (channel_2a.imag.astype(np.int16)).tofile(raw_2a_q)
    with open(os.path.join(root_directory, time_stamp + '_2b_i_.raw'), 'wb') as raw_2b_i:
      (channel_2b.real.astype(np.int16)).tofile(raw_2b_i)
    with open(os.path.join(root_directory, time_stamp + '_2b_q_.raw'), 'wb') as raw_2b_q:
      (channel_2b.imag.astype(np.int16)).tofile(raw_2b_q)

  elif (switch_mode == 3):
    # four polarizations recorded
    channel_1a = data[:,:,0]
    channel_1b = data[:,:,1]
    channel_2a = data[:,:,2]
    channel_2b = data[:,:,3]

    with open(os.path.join(root_directory, time_stamp + '_1a_i_.raw'), 'wb') as raw_1a_i:
      (channel_1a.real.astype(np.int16)).tofile(raw_1a_i)
    with open(os.path.join(root_directory, time_stamp + '_1a_q_.raw'), 'wb') as raw_1a_q:
      (channel_1a.imag.astype(np.int16)).tofile(raw_1a_q)
    with open(os.path.join(root_directory, time_stamp + '_1b_i_.raw'), 'wb') as raw_1b_i:
      (channel_1b.real.astype(np.int16)).tofile(raw_1b_i)
    with open(os.path.join(root_directory, time_stamp + '_1b_q_.raw'), 'wb') as raw_1b_q:
      (channel_1b.imag.astype(np.int16)).tofile(raw_1b_q)

    with open(os.path.join(root_directory, time_stamp + '_2a_i_.raw'), 'wb') as raw_2a_i:
      (channel_2a.real.astype(np.int16)).tofile(raw_2a_i)
    with open(os.path.join(root_directory, time_stamp + '_2a_q_.raw'), 'wb') as raw_2a_q:
      (channel_2a.imag.astype(np.int16)).tofile(raw_2a_q)
    with open(os.path.join(root_directory, time_stamp + '_2b_i_.raw'), 'wb') as raw_2b_i:
      (channel_2b.real.astype(np.int16)).tofile(raw_2b_i)
    with open(os.path.join(root_directory, time_stamp + '_2b_q_.raw'), 'wb') as raw_2b_q:
      (channel_2b.imag.astype(np.int16)).tofile(raw_2b_q)


def range_doppler(data, root_directory, time_stamp, switch_mode, n_seconds, sample_range, presummed_prf, max_range):
  '''
  Produces and saves a Range_Doppler plot of SAR data
  Convert fast-time (axis=0) to range and slow-time (axis=1) to Doppler by way
  of the Fast Fourier Transform
  '''

  #TODO: compute the appropriate axes to display range against Doppler frequency
  # x_axis = np.linspace(0, n_seconds, (int)(n_seconds*presummed_prf) , endpoint=False)  # Doppler axis
  # y_axis = np.linspace(0, max_range, sample_range, endpoint=False)  # range axis

  # select which channel to plot depending on switch mode
  if switch_mode == 1:
    # unpack data into the relevent polarizations depending on hardware configuration
    channel_1a = data[:,:,0]
    channel_1b = data[:,:,1]

    # FFT in the fast-time axis
    channel_1a_fft = np.fft.fft(channel_1a, n=2048, axis=0)
    # FFt in the slow-time axis
    channel_1a_rd = np.fft.fft(np.fft.fft(channel_1a_fft, n=1024, axis=1))

    # FFT in the fast-time axis
    channel_1b_fft = np.fft.fft(channel_1b, n=2048, axis=0)
    # FFt in the slow-time axis
    channel_1b_rd = np.fft.fft(np.fft.fft(channel_1b_fft, n=1024, axis=1))

    # produce the range-Doppler plots
    plot_rd_map(title='rd', data=channel_1a_rd, switch_mode=switch_mode, \
                channel='a', root_directory=root_directory, time_stamp=time_stamp)
    
    plot_rd_map(title='rd', data=channel_1b_rd, switch_mode=switch_mode, \
                channel='b', root_directory=root_directory, time_stamp=time_stamp)
  
  elif switch_mode == 2:
    # unpack data into the relevent polarizations depending on hardware configuration
    channel_2a = data[:,:,0]
    channel_2b = data[:,:,1]

    # FFT in the fast-time axis
    channel_2a_fft = np.fft.fft(channel_2a, n=2048, axis=0)
    # FFt in the slow-time axis
    channel_2a_rd = np.fft.fft(np.fft.fft(channel_2a_fft, n=1024, axis=1))

    # FFT in the fast-time axis
    channel_2b_fft = np.fft.fft(channel_2b, n=2048, axis=0)
    # FFt in the slow-time axis
    channel_2b_rd = np.fft.fft(np.fft.fft(channel_2b_fft, n=1024, axis=1))

    # produce the range-Doppler plots
    plot_rd_map(title='rd', data=channel_2a_rd, switch_mode=switch_mode, \
                channel='a', root_directory=root_directory, time_stamp=time_stamp)
    
    plot_rd_map(title='rd', data=channel_2b_rd, switch_mode=switch_mode, \
                channel='b', root_directory=root_directory, time_stamp=time_stamp)

  elif switch_mode == 3:
    channel_1a = data[:,:,0]
    channel_1b = data[:,:,1]
    channel_2a = data[:,:,2]
    channel_2b = data[:,:,3]

    # FFT in the fast-time axis
    channel_1a_fft = np.fft.fft(channel_1a, n=2048, axis=0)
    # FFt in the slow-time axis
    channel_1a_rd = np.fft.fft(np.fft.fft(channel_1a_fft, n=1024, axis=1))

    # FFT in the fast-time axis
    channel_1b_fft = np.fft.fft(channel_1b, n=2048, axis=0)
    # FFt in the slow-time axis
    channel_1b_rd = np.fft.fft(np.fft.fft(channel_1b_fft, n=1024, axis=1))

    # FFT in the fast-time axis
    channel_2a_fft = np.fft.fft(channel_2a, n=2048, axis=0)
    # FFt in the slow-time axis
    channel_2a_rd = np.fft.fft(np.fft.fft(channel_2a_fft, n=1024, axis=1))

    # FFT in the fast-time axis
    channel_2b_fft = np.fft.fft(channel_2b, n=2048, axis=0)
    # FFt in the slow-time axis
    channel_2b_rd = np.fft.fft(np.fft.fft(channel_2b_fft, n=1024, axis=1))

    # produce the range-Doppler plots for each polarization
    plot_rd_map(title='rd', data=channel_1a_rd, switch_mode=1, \
                channel='a', root_directory=root_directory, time_stamp=time_stamp)
    
    plot_rd_map(title='rd', data=channel_1b_rd, switch_mode=1, \
                channel='b', root_directory=root_directory, time_stamp=time_stamp)
          
    plot_rd_map(title='rd', data=channel_2a_rd, switch_mode=2, \
                channel='a', root_directory=root_directory, time_stamp=time_stamp)
    
    plot_rd_map(title='rd', data=channel_2b_rd, switch_mode=2, \
                channel='b', root_directory=root_directory, time_stamp=time_stamp)


def plot_rd_map(title, data, switch_mode, channel, root_directory, time_stamp):
  '''
  Plot Range-Doppler map figure to file
  '''

  plt.subplots(nrows=1,ncols=1)
  plt.title(title.upper() + " " + str(switch_mode) + channel)
  plt.xlabel("Doppler Bins")
  plt.ylabel("Range Bins")
    
  plt.imshow(20*np.log10(np.abs(data)), aspect='auto', origin='lower')
    
  file_name = time_stamp + "_" + title + "_" + str(switch_mode) + channel + ".png"
  plt.savefig(os.path.join(root_directory, 'quicklook/'+file_name))
  
  print(time_stamp + " " + title.upper() + " " + str(switch_mode) + channel + " saved.")


def rti(data, root_directory, time_stamp, switch_mode, n_seconds, n_pri, presummed_prf, max_range):
  '''
  Produces and saves a range time intensity plot of the data
  Compute a 2D Fast Fourier Transform to convert fast-time axis to range and slow-time axis to time
  '''
  #TODO: Reduce dynamic range of image in order to see more low intensity clutter returns

  # Compute the axes
  x_axis = np.linspace(0, n_seconds, n_pri, endpoint=False)  # time axis
  y_axis = np.linspace(0, max_range, (int)(n_seconds*presummed_prf), endpoint=False)  # range axis

  # select which channel to plot depending on switch mode
  if switch_mode == 1:
    channel_1a = data[:,:,0]
    channel_1b = data[:,:,1]
    channel_1a_fft = np.fft.fftshift(np.fft.fft2(channel_1a))
    channel_1b_fft = np.fft.fftshift(np.fft.fft2(channel_1b))

    plot_rti(title='rti', data=channel_1a_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=switch_mode, \
              channel='a', root_directory=root_directory, time_stamp=time_stamp)
    
    plot_rti(title='rti', data=channel_1b_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=switch_mode, \
              channel='b', root_directory=root_directory, time_stamp=time_stamp)
  
  elif switch_mode == 2:
    channel_2a = data[:,:,0]
    channel_2b = data[:,:,1]
    channel_2a_fft = np.fft.fftshift(np.fft.fft2(channel_2a))
    channel_2b_fft = np.fft.fftshift(np.fft.fft2(channel_2b))

    plot_rti(title='rti', data=channel_2a_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=switch_mode, \
              channel='a', root_directory=root_directory, time_stamp=time_stamp)
    
    plot_rti(title='rti', data=channel_2b_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=switch_mode, \
              channel='b', root_directory=root_directory, time_stamp=time_stamp)

  elif switch_mode == 3:
    channel_1a = data[:,:,0]
    channel_1b = data[:,:,1]
    channel_2a = data[:,:,2]
    channel_2b = data[:,:,3]
    
    channel_1a_fft = np.fft.fftshift(np.fft.fft2(channel_1a))
    channel_1b_fft = np.fft.fftshift(np.fft.fft2(channel_1b))
    channel_2a_fft = np.fft.fftshift(np.fft.fft2(channel_2a))
    channel_2b_fft = np.fft.fftshift(np.fft.fft2(channel_2b))

    plot_rti(title='rti', data=channel_1a_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=1, \
              channel='a', root_directory=root_directory, time_stamp=time_stamp)
    
    plot_rti(title='rti', data=channel_1b_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=1, \
              channel='b', root_directory=root_directory, time_stamp=time_stamp)
          
    plot_rti(title='rti', data=channel_2a_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=2, \
              channel='a', root_directory=root_directory, time_stamp=time_stamp)
    
    plot_rti(title='rti', data=channel_2b_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=2, \
              channel='b', root_directory=root_directory, time_stamp=time_stamp)


def plot_rti(title, data, x_axis, y_axis, switch_mode, channel, root_directory, time_stamp):
  '''
  Plot and save Range-Time Intensity to file
  '''

  plt.subplots(nrows=1,ncols=1)
  plt.title(title.upper() + " " + str(switch_mode) + channel)
  plt.xlabel("Slow Time [s]")
  plt.ylabel("Range [m]")
  
  data_log10 = 20*np.log10(np.abs(data)) - np.amax(20*np.log10(np.abs(data)))
  
  plt.imshow(data_log10/1.1, aspect='auto', origin='lower', extent=[min(x_axis), max(x_axis), min(y_axis), max(y_axis)])
    
  file_name = time_stamp + "_" + title + "_" + str(switch_mode) + channel + ".png"
  plt.savefig(os.path.join(root_directory, 'quicklook/'+file_name))
  
  print(time_stamp + " " + title.upper() + " " + str(switch_mode) + channel + " saved.")


def plot_histogram(data, root_directory):
  '''
  Plot a histogram of the data to visualise the spread of data
  Not currently working as intended
  '''
  #TODO: Determine what information to plot that will make a meaning report of the data

  channel_1a = data[:,:,0]
  channel_1b = data[:,:,1]
  channel_1a_fft = np.fft.fftshift(np.fft.fft(channel_1a, axis=0))
  channel_1b_fft = np.fft.fftshift(np.fft.fft(channel_1b, axis=0))
  
  plt.figure()
  plt.hist(np.abs(channel_1a_fft).flatten(), bins='auto')
  plt.savefig(os.path.join(root_directory, "data_histogram_abs_1a.png"))

  plt.figure()
  plt.hist(np.abs(channel_1b_fft).flatten(), bins='auto')
  plt.savefig(os.path.join(root_directory, "data_histogram_abs_1b.png"))

  print("Histograms saved.")


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
  if dlat > 10e-9:    # provide an error margin
    distance = -1*R*c*1000
  if dlat < -10e-9:
    distance = R*c*1000

  return distance # in metres


def motion_compensation(data, x_axis, switch_mode, root_directory):
  '''
  Read motion data from JSON generated by system
  Extract the following information:
    - timestamps, timestamp
    - latitude and longitude, lat & long
    - altitude, alt
    - absolute velocity, abs_vel
    - expected aperture, e
    - offset normal, n
  '''
  motion_filename = glob.glob(os.path.join(root_directory, '*.json'))[0]
  
  # read motion file contents
  motion_data = [json.loads(line) for line in open(motion_filename, 'r')]

  timestamp = np.array([])
  altitude = np.array([])
  velocity = np.array([])
  aperture = np.array([])
  normal = np.array([])
  latitude = np.array([])
  longitude = np.array([])

  for line in motion_data:
    timestamp = np.append(timestamp, line['timestamp'])
    altitude = np.append(altitude, line['u'])
    velocity = np.append(velocity, line['abs_vel'])
    aperture = np.append(aperture, line['e'])
    normal = np.append(normal, line['n'])
    latitude = np.append(latitude, line['lat'])
    longitude = np.append(longitude, line['lon'])
  
  # form nominal flight path from start and end lat-long coordinates - straight line fit
  fit = np.polyfit(np.array([longitude[0], longitude[-1]]), np.array([latitude[0], latitude[-1]]), deg=1)
  
  # nominal path will be a line through (longitude, nom_path(longitude))
  nom_path = np.poly1d(fit) # the nominal path the plane should have 
  
  ###### test statement ######
  plt.figure()
  plt.plot(longitude, latitude, '-', label='Actual')
  plt.plot(longitude, nom_path(longitude), '--', label='Nominal')
  plt.xlabel('Longitude')
  plt.ylabel('Latitude')
  plt.title('Nominal vs Actual Flight Path')
  plt.legend()
  plt.grid()
  plt.savefig(os.path.join(root_directory, 'quicklook/actual_vs_nomimal_path.png'))

  # now to find the deviations of the flight path from the straight-line simulated above
  lon_A = longitude
  lon_B = longitude
  lat_A = nom_path(longitude)
  lat_B = latitude
  range_deviation = np.array([])
  for i in range(longitude.size):
    range_deviation  = np.append(range_deviation, haversine(lon_A[i], lat_A[i], lon_B[i], lat_B[i]))

  # resize range deviation map to number of range bins
  n_pris = data.shape[0]
  n_range_bins = data.shape[1]
  n_polarizations = data.shape[2]

  t_interpolate = np.array([])
  for t in timestamp:
    t_interpolate = np.append(t_interpolate, (t-timestamp[0])/1e3)  # given 10Hz refresh rate, t has 0.1s increments

  range_dev = interpolate.interp1d(t_interpolate, range_deviation, kind='linear')
  range_deviation_correction = range_dev(x_axis)
  
  # apply range bin correction to the data
  # loop through each range bin in each pri and apply range correction
  phase_shift = np.exp(np.multiply(-1j*(4*np.pi), range_deviation_correction))
  phase_shift = np.array([phase_shift,])
  data_out = np.array([])

  if switch_mode  == 1:
    channel_1a = data[:,:,0]
    channel_1b = data[:,:,1]

    temp_1a = np.multiply(np.transpose(channel_1a), phase_shift)
    temp_1b = np.multiply(np.transpose(channel_1b), phase_shift)
    
    data_out = np.dstack((temp_1a, temp_1b))
  
  if switch_mode  == 2:
    channel_2a = data[:,:,0]
    channel_2b = data[:,:,1]

    temp_2a = np.multiply(np.transpose(channel_2a), phase_shift)
    temp_2b = np.multiply(np.transpose(channel_2b), phase_shift)

    # TODO: Range bin shifting
    # apply range bin shifting to the dataset
    #for pri in range(n_pris):
      #temp_2a[pri,:] = temp_2a[pri,:] + range_deviation_correction
      #temp_2b[pri,:] = temp_2b[pri,:] + range_deviation_correction

    data_out = np.dstack((temp_2a, temp_2b))
        
  elif switch_mode == 3:
    channel_1a = data[:,:,0]
    channel_1b = data[:,:,1]
    channel_2a = data[:,:,2]
    channel_2b = data[:,:,3]

    temp_1a = np.multiply(np.transpose(channel_1a), phase_shift)
    temp_1b = np.multiply(np.transpose(channel_1b), phase_shift)
    temp_2a = np.multiply(np.transpose(channel_2a), phase_shift)
    temp_2b = np.multiply(np.transpose(channel_2b), phase_shift)
 
    data_out = np.dstack((temp_1a,    # 1a
                          temp_1b,    # 1b
                          temp_2a,    # 2a
                          temp_2b))   # 2b

  return data_out


if __name__ == "__main__":
  main()  # call the main function