"""
This is a quicklook processor which is used to give a quick
snapshot of miloSAR data for review

Author:   Tlotliso Mapana
Date:     28 September 2021
"""

from cProfile import label
from cmath import log10, phase
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
from string import *
from PIL import Image

import g2tools
#from g2main import *
#from g2tools import *

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
    time_stamp, switch_mode, n_seconds, presummed_prf, n_pri, min_range, max_range, \
          n_range_bins, sampling_rate, n_az_points \
          = load_configuration_parameters(root_directory=root_directory)

    # load data from file and pack into a data matrix
    data = load_data(root_directory=root_directory, time_stamp=time_stamp, switch_mode=switch_mode, \
                    n_seconds=n_seconds, presummed_prf=presummed_prf, n_pri=n_pri)

    # motion compensation
    x_axis = np.linspace(0, n_seconds, n_pri, endpoint=False)  # time axis
    y_axis = np.linspace(0, max_range, (int)(n_seconds*presummed_prf), endpoint=False)  # range axis

    # OPTIONS
    if len(sys.argv) > 2:
      for opt in range(2, len(sys.argv)):
        if sys.argv[opt] == 'raw':
          # save binary image to file
          save_raw(data, root_directory, switch_mode, n_range_bins, n_seconds, time_stamp)

        elif sys.argv[opt] == 'rd':
          # produce a range-Doppler plot and save to directory
          range_doppler(data, root_directory, time_stamp, switch_mode, \
                        n_seconds, presummed_prf, max_range, sampling_rate, all=False)

        elif sys.argv[opt] == 'rti':
          # produce range-time plot and save to diectory
          rti(data, root_directory, time_stamp, switch_mode, \
              n_seconds, n_pri, presummed_prf, max_range, all=False)
    
        elif sys.argv[opt] == 'hist':
          # produce a histogram of the raw data
          plot_histogram(data, root_directory=root_directory)

        elif sys.argv[opt] == 'mocomp':
          # perform motion compensation
          data = motion_compensation(data=data, x_axis=x_axis, switch_mode=switch_mode, \
                                      root_directory=root_directory, plot_path=True)
        
        elif sys.argv[opt] == 'sar':
          # run SAR processor
          data = fast_time_fft(data=data)          
          SAR(data=data, root_directory=root_directory, switch_mode=switch_mode, \
              time_stamp=time_stamp, prf=presummed_prf, n_az_points=n_az_points, \
              n_range_bins=n_range_bins, max_range=max_range, n_seconds=n_seconds, \
              d_range=[])

        else:
          print("\nUnidentified function selected.")
          print("Usage: processor.py dataset_root_directory [OPTIONS]")
          print("  OPTIONS:")
          print("    raw:    Save image in raw binary file")
          print("    rd:     Plot range-Doppler map for each polarization")
          print("    rti:    Plot and save range time intensity plot for each polarization")
          print("    sar:    Compute and save SAR image\n")

    else:
      print("\nUsage: processor.py dataset_root_directory [OPTIONS]")
      print("  OPTIONS:")
      print("    raw:    Save image in raw binary file")
      print("    rd:     Plot range-Doppler map for each polarization")
      print("    rti:    Plot and save range time intensity plot for each polarization")
      print("    sar:    Compute and save SAR image\n")
              
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

  min_range_string = setup['geometry']['min_range']  # 100.0; [m]
  min_range = int(min_range_string.split('.')[0])
  max_range_string = setup['geometry']['max_range']  # 900.0; [m]
  max_range = int(max_range_string.split('.')[0])

  n_range_bins = max_range - min_range

  # realisable prf of data file stored in memory as data rate is halved when writing to memory
  presummed_prf = prf/presumming_factor # TODO rename to more understandable term
  #if switch_mode == 3:
    #presummed_prf = presummed_prf/2

  n_az_points = int(n_seconds*presummed_prf)
  if switch_mode == 3:
    n_az_points = int(n_seconds*presummed_prf//2)

  # number of range bins recorded for each pri
  n_pri = end_index - start_index + 1
  
  print("Configuration parameters loaded.")

  return time_stamp, switch_mode, n_seconds, presummed_prf, n_pri, min_range, max_range, n_range_bins, sampling_rate, n_az_points


def load_data(root_directory, time_stamp, switch_mode, n_seconds, presummed_prf, n_pri):
  '''
  This function reads in raw miloSAR data stored in a .bin file in the given directory
  '''
  # open file stored at existing directory and read as binary
  file_name = glob.glob(os.path.join(root_directory, '*.bin'))[0]  
  with open(file_name, 'rb') as f:
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
  n_az_points = int(n_seconds*presummed_prf)
  n_samples_per_channel = int(n_az_points*n_pri)

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
    n_az_points, n_pri))
  channel_b_data = np.transpose(channel_b_data.reshape(
    n_az_points, n_pri))
  
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

  #TODO: Constants
  rp_adc_vpp = 2 # 2V peak-to-peak on this model red pitaya
  rp_adc_bits = 14 # for this model
  adc_bits_volts = 2*rp_adc_vpp/math.pow(2, rp_adc_bits)

  # scale data
  data = data*adc_bits_volts

  # apply windowing function to data
  data = window(data=data, switch_mode=switch_mode, n_pri=n_pri)
    
  print("Dataset", time_stamp, "loaded successfully.")
  return data


def window(data, switch_mode, n_pri):
    '''
    Attempt at windowing the data to remove high sidelobes
    Using Hamming window function
    '''
    window_function = np.hamming(n_pri)

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


def fast_time_fft(data):
  n_az_points = data.shape[0]
  n_FFT = int(math.pow(2, math.ceil(np.log2(n_az_points))))
  data = np.fft.fftshift(np.fft.fft(data, n_FFT, axis=0), axes=0)

  return data


def save_raw(data, root_directory, switch_mode, n_range_bins, n_seconds, time_stamp):
  '''
  Not currently saving raw data as intended
  TODO: Rewrite this function to produce correct results
  '''
  data = fast_time_fft(data)
  image = pow(np.mean(abs(data), axis=2), 2)
  image.astype('complex64').tofile(
    os.path.join(root_directory, 'quicklook/' + time_stamp + '.bin') )

  data_dB = 20*np.log10(np.abs(image))

  plt.subplots(nrows=1,ncols=1)
  plt.xlabel("Azimuth")
    
  plt.imshow(np.abs(data_dB), cmap='gray', aspect='equal', origin='upper', vmin=None, vmax=None)
  img_file = time_stamp + "_" + str(switch_mode) + ".png"
  image_path = os.path.join(root_directory, 'quicklook/'+img_file)
  plt.imsave(image_path, np.abs(data_dB), cmap='gray', origin='upper', format='png')

  synthetic_aperture = n_seconds*30

  # pixel scaling for SAR image
  scaled_width = int(n_range_bins / 900 * synthetic_aperture)
  scaled_height = n_range_bins

  image = Image.open(image_path)
  image.load()
  image = image.resize((scaled_width, scaled_height), resample=Image.LANCZOS)
  image.save(os.path.splitext(image_path)[
            0] + os.path.splitext(image_path)[1], format='png')

  print('Binary image save to ' + image_path)


def range_doppler(data, root_directory, time_stamp, switch_mode, n_seconds, presummed_prf, max_range, sampling_rate, all=False):
  '''
  Produces and saves a Range_Doppler plot of SAR data
  Convert fast-time (axis=1) to range and slow-time (axis=0) to Doppler by way
  of the Fast Fourier Transform
  '''

  # Set FFT length
  n_FFT_Doppler = 64

  #TODO: compute the appropriate axes to display range against Doppler frequency
  x_axis = np.fft.fftfreq(n=n_FFT_Doppler, d=1/sampling_rate)/1e9 # frequecy axis scaled
  y_axis = np.linspace(0, max_range, (int)(n_seconds*presummed_prf), endpoint=False)  # range axis

  # select which channel to plot depending on switch mode
  if switch_mode == 1:
    # unpack data into the relevent polarizations depending on hardware configuration
    channel_1a = data[:,:,0]
    channel_1b = data[:,:,1]

    # FFT in the slow-time axis
    channel_1a_fft = np.fft.fft(channel_1a, n=n_FFT_Doppler, axis=1)
    # produce the range-Doppler plot
    plot_rd_map(title='rd', data=channel_1a_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=switch_mode, \
                channel='a', root_directory=root_directory, time_stamp=time_stamp)

    if all:
      # FFT in the slow-time axis
      channel_1b_fft = np.fft.fft(channel_1b, n=n_FFT_Doppler, axis=1)
      # produce the range-Doppler plot
      plot_rd_map(title='rd', data=channel_1b_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=switch_mode, \
                channel='b', root_directory=root_directory, time_stamp=time_stamp)
  
  elif switch_mode == 2:
    # unpack data into the relevent polarizations depending on hardware configuration
    channel_2a = data[:,:,0]
    channel_2b = data[:,:,1]

    # FFT in the slow-time axis
    channel_2a_fft = np.fft.fft(channel_2a, n=n_FFT_Doppler, axis=1)
    # produce the range-Doppler plot
    plot_rd_map(title='rd', data=channel_2a_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=switch_mode, \
                channel='a', root_directory=root_directory, time_stamp=time_stamp)
    
    if all:
      # FFT in the slow-time axis
      channel_2b_fft = np.fft.fft(channel_2b, n=n_FFT_Doppler, axis=1)
      # produce the range-Doppler plot
      plot_rd_map(title='rd', data=channel_2b_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=switch_mode, \
                channel='b', root_directory=root_directory, time_stamp=time_stamp)

  elif switch_mode == 3:
    channel_1a = data[:,:,0]
    channel_1b = data[:,:,1]
    channel_2a = data[:,:,2]
    channel_2b = data[:,:,3]

    # FFT in the slow-time axis
    channel_1a_fft = np.fft.fft(channel_1a, n=n_FFT_Doppler, axis=1)
    # produce the range-Doppler plot
    plot_rd_map(title='rd', data=channel_1a_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=1, \
                channel='a', root_directory=root_directory, time_stamp=time_stamp)

    if all:
      # FFT in the slow-time axis
      channel_1b_fft = np.fft.fft(channel_1b, n=n_FFT_Doppler, axis=1)
      channel_2a_fft = np.fft.fft(channel_2a, n=n_FFT_Doppler, axis=1)
      channel_2b_fft = np.fft.fft(channel_2b, n=n_FFT_Doppler, axis=1)
      
      # produce the range-Doppler plot
    
      plot_rd_map(title='rd', data=channel_1b_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=1, \
                channel='b', root_directory=root_directory, time_stamp=time_stamp)
          
      plot_rd_map(title='rd', data=channel_2a_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=2, \
                channel='a', root_directory=root_directory, time_stamp=time_stamp)
    
      plot_rd_map(title='rd', data=channel_2b_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=2, \
                channel='b', root_directory=root_directory, time_stamp=time_stamp)


def plot_rd_map(title, data, x_axis, y_axis, switch_mode, channel, root_directory, time_stamp):
  '''
  Plot Range-Doppler map figure to file
  '''

  plt.subplots(nrows=1,ncols=1)
  plt.title(title.upper() + " " + str(switch_mode) + channel)
  plt.xlabel("Frequency [GHz]")
  plt.ylabel("Range [m]")

  data_dB = 20*np.log10(np.abs(data))
  data_max = np.amax(data_dB)
  data_min = data_max - 30
  data_dB = np.clip(data_dB, data_min, data_max)
    
  plt.imshow(X=data_dB, aspect='auto', origin='lower', extent=[min(x_axis), max(x_axis), min(y_axis), max(y_axis)])
  plt.colorbar(label='Power [dB]')
  plt.tight_layout()
    
  file_name = time_stamp + "_" + title + "_" + str(switch_mode) + channel + ".png"
  plt.savefig(os.path.join(root_directory, 'quicklook/'+file_name))
  
  print(time_stamp + " " + title.upper() + " " + str(switch_mode) + channel + " saved.")


def rti(data, root_directory, time_stamp, switch_mode, n_seconds, n_pri, presummed_prf, max_range, all=False):
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
    # unpack data into the relevent polarizations depending on hardware configuration
    channel_1a = data[:,:,0]
    channel_1b = data[:,:,1]

    channel_1a_fft = fast_time_fft(data=channel_1a)
    plot_rti(title='rti', data=channel_1a_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=switch_mode, \
              channel='a', root_directory=root_directory, time_stamp=time_stamp)

    if all:
      channel_1b_fft = fast_time_fft(data=channel_1b)
      plot_rti(title='rti', data=channel_1b_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=switch_mode, \
              channel='b', root_directory=root_directory, time_stamp=time_stamp)
  
  elif switch_mode == 2:
    # unpack data into the relevent polarizations depending on hardware configuration
    channel_2a = data[:,:,0]
    channel_2b = data[:,:,1]

    # FFT in the fast-time axis
    channel_2a_fft = fast_time_fft(data=channel_2a)
    plot_rti(title='rti', data=channel_2a_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=switch_mode, \
              channel='a', root_directory=root_directory, time_stamp=time_stamp)

    if all:
      # FFT in the fast-time axis
      channel_2b_fft = fast_time_fft(data=channel_2b)
      plot_rti(title='rti', data=channel_2b_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=switch_mode, \
              channel='b', root_directory=root_directory, time_stamp=time_stamp)

  elif switch_mode == 3:
    channel_1a = data[:,:,0]
    channel_1b = data[:,:,1]
    channel_2a = data[:,:,2]
    channel_2b = data[:,:,3]
    
    channel_1a_fft = fast_time_fft(data=channel_1a)
    plot_rti(title='rti', data=channel_1a_fft, x_axis=x_axis, y_axis=y_axis, switch_mode=1, \
              channel='a', root_directory=root_directory, time_stamp=time_stamp)
    
    if all:
      channel_1b_fft = fast_time_fft(data=channel_1b)
      channel_2a_fft = fast_time_fft(data=channel_2a)
      channel_2b_fft = fast_time_fft(data=channel_2b)
      
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
  
  data_dB = 20*np.log10(abs(data))
  data_max = np.amax(data_dB)
  data_min = np.nanmin(np.mean(data_dB, axis=1))

  data_dB = np.clip(data_dB, data_min, data_max)
  
  plt.imshow(data_dB, interpolation='none', cmap='viridis', aspect='auto', origin='lower', extent=[min(x_axis), max(x_axis), min(y_axis), max(y_axis)])
    
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
  if dlat > 0:    # provide an error margin
    distance = -1*R*c*1000
  if dlat < 0:
    distance = 1*R*c*1000

  return distance # in metres


def motion_compensation(data, x_axis, switch_mode, root_directory, plot_path=False):
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
  
  if plot_path:
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

  t_interpolate = np.array([])
  for t in timestamp:
    t_interpolate = np.append(t_interpolate, (t-timestamp[0])/1e3)  # given 10Hz refresh rate, t has 0.1s increments

  range_dev = interpolate.interp1d(t_interpolate, range_deviation, kind='linear')
  x_axis = np.linspace(0, 30, 2048, endpoint=False)  # time axis
  range_deviation_correction = range_dev(x_axis)
  
  # apply range bin correction to the data
  # loop through each range bin in each pri and apply range correction
  phase_shift = np.exp(np.multiply(-1j*(4*np.pi)/2437498854.473165, range_deviation_correction))
  phase_shift = np.array([phase_shift,])
  phase_shift = np.transpose(phase_shift)

  data_fft = np.fft.fft(data, 2048, axis=0)
  data_out = np.array([])

  if switch_mode  == 1:
    channel_1a = data_fft[:,:,0]
    channel_1b = data_fft[:,:,1]

    temp_1a = np.multiply(channel_1a, phase_shift)
    temp_1a = np.fft.ifft(temp_1a, axis=0)

    temp_1b = np.multiply(channel_1b, phase_shift)
    temp_1b = np.fft.ifft(temp_1b, axis=0)
    
    data_out = np.dstack((temp_1a, temp_1b))
  
  if switch_mode  == 2:
    channel_2a = data_fft[:,:,0]
    channel_2b = data_fft[:,:,1]
    
    temp_2a = np.multiply(channel_2a, phase_shift)
    temp_2a = np.fft.ifft(temp_2a, axis=0)

    temp_2b = np.multiply(channel_2b, phase_shift)
    temp_2b = np.fft.ifft(temp_2b, axis=0)

    data_out = np.dstack((temp_2a, temp_2b))
        
  elif switch_mode == 3:
    channel_1a = data_fft[:,:,0]
    channel_1b = data_fft[:,:,1]
    channel_2a = data_fft[:,:,2]
    channel_2b = data_fft[:,:,3]

    temp_1a = np.multiply(channel_1a, phase_shift)
    temp_1a = np.fft.ifft(temp_1a, axis=0)

    temp_1b = np.multiply(channel_1b, phase_shift)
    temp_1b = np.fft.ifft(temp_1b, axis=0)

    temp_2a = np.multiply(channel_2a, phase_shift)
    temp_2a = np.fft.ifft(temp_2a, axis=0)

    temp_2b = np.multiply(channel_2b, phase_shift)
    temp_2b = np.fft.ifft(temp_2b, axis=0)
 
    data_out = np.dstack((temp_1a,    # 1a
                          temp_1b,    # 1b
                          temp_2a,    # 2a
                          temp_2b))   # 2b
  
  #print(data_out.shape) # debug
  
  return data_out


def SAR(data, root_directory, switch_mode, time_stamp, prf, n_az_points, n_range_bins, max_range, n_seconds, d_range=[]):
  '''
  Make call to G2 program
  '''
  # TODO: for now
  synthetic_aperture = n_seconds*30

  # pixel scaling for SAR image
  scaled_width = int(n_range_bins / max_range * synthetic_aperture)
  scaled_height = n_range_bins
  G2_out = data

  if switch_mode == 1:
    G2_out_1a = G2(data=data[:,:,0], n_az_points=n_az_points, title='1a', root_directory=root_directory, time_stamp=time_stamp, prf=prf, n_range_bins=n_range_bins)
    G2_out_1b = G2(data=data[:,:,1], n_az_points=n_az_points, title='1b', root_directory=root_directory, time_stamp=time_stamp, prf=prf, n_range_bins=n_range_bins)

    G2_out = np.dstack((G2_out_1a, G2_out_1b))

  elif switch_mode == 2:
    G2_out_2a = G2(data=data[:,:,0], n_az_points=n_az_points, title='2a', root_directory=root_directory, time_stamp=time_stamp, prf=prf, n_range_bins=n_range_bins)
    G2_out_2b = G2(data=data[:,:,1], n_az_points=n_az_points, title='2b', root_directory=root_directory, time_stamp=time_stamp, prf=prf, n_range_bins=n_range_bins)

    G2_out = np.dstack((G2_out_2a, G2_out_2b))
  
  elif switch_mode == 3:
    
    G2_out_1a = G2(data=data[:,:,0], n_az_points=n_az_points, title='1a', root_directory=root_directory, time_stamp=time_stamp, prf=prf, n_range_bins=2000)
    G2_out_1b = G2(data=data[:,:,1], n_az_points=n_az_points, title='1b', root_directory=root_directory, time_stamp=time_stamp, prf=prf, n_range_bins=2000)
    G2_out_2a = G2(data=data[:,:,2], n_az_points=n_az_points, title='2a', root_directory=root_directory, time_stamp=time_stamp, prf=prf, n_range_bins=2000)
    G2_out_2b = G2(data=data[:,:,3], n_az_points=n_az_points, title='2b', root_directory=root_directory, time_stamp=time_stamp, prf=prf, n_range_bins=2000)

    G2_out= np.dstack((G2_out_1a, G2_out_1b))
    G2_out = np.dstack((G2_out, G2_out_2a))
    G2_out = np.dstack((G2_out, G2_out_2b))
  
  # noncoherent avering of data
  image = pow(np.mean(abs(G2_out), axis=2), 2)
  image.astype('complex64').tofile(
			os.path.join(root_directory, 'quicklook/' + 'image.bin') )
    
  plt.subplots(nrows=1,ncols=1)
  plt.xlabel('Azimuth [m]')
  plt.ylabel('Range [m]')

  image = 20*np.log10(abs(image))
  img_mean = np.mean(image, axis=1)
  
  # dynamic range
  a_min = np.nanmin(img_mean[img_mean!=-np.inf]) + 50
  a_max = np.amax(image)
  if d_range:
    a_min = d_range[0]
    a_min = d_range[1]
  
  print(a_max)
  image = np.clip(image, a_min=a_min, a_max=a_max)
    
  plt.imshow(image, cmap='gray', aspect='equal', origin='upper', vmin=None, vmax=None)
  file_name = time_stamp + '.sar.noco'
  image_path = os.path.join(root_directory, 'quicklook/'+file_name)
  plt.imsave(image_path, image, cmap='gray', origin='upper', format='png')

  sar_image = Image.open(image_path)
  sar_image.load()
  sar_image = sar_image.resize((scaled_width, scaled_height), resample=Image.LANCZOS)
  sar_image.save(os.path.splitext(image_path)[
               0] + os.path.splitext(image_path)[1], format='png')
    
  print("SAR image saved.")


def G2(data, n_az_points, title, root_directory, time_stamp, prf, n_range_bins):
  '''
  SAR processor using G2
  (C) J Horrell (1999)
  '''
  
  # create g2 files required
  g2_input = os.path.join(root_directory, 'quicklook/' + \
          time_stamp + '_' + title + '.rng')
  g2_cmd = os.path.join(root_directory, 'quicklook/' + \
          time_stamp + '_' + title + '.cmd')
  g2_log = os.path.join(root_directory, 'quicklook/' + \
          time_stamp + '_' + title + '.log')
  g2_out = os.path.join(root_directory, 'quicklook/' + \
          time_stamp + '_' + title + '.azi')
  
  # save data to file
  data.astype('complex64').tofile(g2_input)

  with open(g2_cmd, 'w') as f_id:
    f_id.write('miloSAR azimuth compression command file (azcom)\n')
    f_id.write('$ProgramVersion (jmh)         => 1.1\n\n')

    f_id.write('$ScreenUpdateRate             => ' + str(1) + '\n')
    f_id.write('$LogFileName                  => ' +
            str(g2_log) + '\n')
    f_id.write('$InputStartSampleDelay        => ' + str(0) + '\n')
    f_id.write('$CarrierFreq [Hz]             => ' + str(2437498854.473165) + '\n')
    f_id.write('$InputPRF [Hz]                => ' + str(prf) + '\n')
    f_id.write('$NomGroundSpeed [m/s]         => ' + str(32) + '\n')
    f_id.write('$InputFileAzPts               => ' + str(n_az_points) + '\n')
    f_id.write('$StartProcessAzPt             => ' + str(0) + '\n') #change
    f_id.write('$AzPtsToProcess               => ' + str(n_az_points) + '\n')
    f_id.write('$InputFileRngBins             => ' + str(n_range_bins) + '\n')
    f_id.write('$StartProcessRngBin           => ' + str(0) + '\n') #change
    f_id.write('$RngBinsToProcess             => ' + str(n_range_bins) + '\n')
    f_id.write('$InputDCOffsetI               => ' + str(0.0) + '\n')
    f_id.write('$InputDCOffsetQ               => ' + str(0.0) + '\n')
    f_id.write('$InvFFTSizeReduc [pow of 2]   => ' + str(1) + '\n')
    f_id.write('$InputFileName                => ' + str(g2_input) + '\n')
    f_id.write('$OutputFileName               => ' + str(g2_out) + '\n')
    f_id.write('$AppendExistOutFileFlg [Y/N]  => ' + str('N') + '\n')
    f_id.write('$RngFocSegments               => ' + str(-1) + '\n')
    f_id.write('$RefFuncSign [+-1]            => ' + str(1) + '\n') #change
    f_id.write('$A2DFreq [Hz]                 => ' + str(284881608.77714235) + '\n')
    f_id.write('$NomAzRes [m]                 => ' + str(0.5475561694721032) + '\n')
    f_id.write('$WinConstTime [0.0-1.0]       => ' + str(0.08) + '\n')
    f_id.write('$NumLooks                     => ' + str(1) + '\n')
    f_id.write('$LookOverlapFrac [0.0-1.0]    => ' + str(0.0) + '\n')
    f_id.write('$WinConstFreq [0.0-1.0]       => ' + str(0.08) + '\n')
    f_id.write('$RngCurvInterpSize            => ' + str(4) + '\n')
    f_id.write('$RngCurvBatchSize             => ' + str(16) + '\n') #play with this
    f_id.write('$PostSumRatio                 => ' + str(1) + '\n')
    f_id.write('$DetectMethod                 => ' + str(0) + '\n')
    f_id.write('$InputDataType                => ' + str(3) + '\n')
    f_id.write('$OutputDataType               => ' + str(3) + '\n')
    f_id.write('$Scale                        => ' + str(1) + '\n')
    f_id.write('$ReportMax [1/0]              => ' + str(1) + '\n')

  azcom_executable = os.path.join(os.getcwd(), 'g2', 'azcom')
  os.system(azcom_executable + ' ' + str(g2_cmd) + ' > ' + root_directory + '/quicklook/g2out.txt')

  image = np.fromfile(g2_out, dtype='complex64')
  #image = np.fromfile(g2_input, dtype='complex64')
  G2_out = np.flipud(image.reshape(n_range_bins, n_az_points))
  #G2_out = np.flipud(image.reshape(2048, n_az_points))
  G2_out = np.nan_to_num(G2_out)

  return G2_out


def range_compression(data, root_directory, time_stamp, plot=True):
  '''
  SAR range compression
  '''
  rnc_data = fast_time_fft(data=data)

  if plot==True:
    plt.subplots(nrows=1,ncols=1)
    plt.title('Azimuth Compresed')
    plt.xlabel("Slow Time [s]")
    plt.ylabel("Range [m]")
    
    data_dB = 20*np.log10(abs(rnc_data))
    data_max = np.amax(data_dB)
    data_min = np.nanmin(np.mean(data_dB, axis=1))

    data_dB = np.clip(data_dB, data_min, data_max)
    
    plt.imshow(data_dB, interpolation='none', cmap='viridis', aspect='auto', origin='lower', extent=[min(x_axis), max(x_axis), min(y_axis), max(y_axis)])
      
    file_name = time_stamp + ".png"
    plt.savefig(os.path.join(root_directory, 'quicklook/'+file_name))


if __name__ == "__main__":
  main()  # call the main function