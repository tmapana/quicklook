"""
This is a quicklook processor which is used to give a quick
snapshot of miloSAR data for review

Author:   Tlotliso Mapana
Date:     28 September 2021
"""

import os
import sys
import math
from matplotlib.colors import LogNorm
import numpy as np
import configparser
from numpy.fft import fftfreq
from numpy.lib.function_base import add_newdoc_ufunc
import matplotlib.pyplot as plt
from numpy.linalg import matrix_rank

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
    time_stamp, switch_mode, n_seconds, presummed_prf, sampling_rate, range_bins, max_range \
       = load_configuration_parameters(root_directory)

    # load data from file and pack into a data matrix
    data = load_data(root_directory, time_stamp, switch_mode, n_seconds, range_bins, presummed_prf)

    # OPTIONS
    if len(sys.argv) > 2:
      for opt in range(2, len(sys.argv)):
        if sys.argv[opt] == 'raw':
          # save .raw file
          save_raw(data, root_directory, time_stamp, switch_mode)
        
        elif sys.argv[opt] == 'rd':
          # produce a range-Doppler plot and save to directory
          range_doppler(data, root_directory, time_stamp, switch_mode, \
                        n_seconds, range_bins, presummed_prf, max_range)

        elif sys.argv[opt] == 'rti':
          # produce range-time plot and save to diectory
          rti(data, root_directory, time_stamp, switch_mode, \
              n_seconds, range_bins, presummed_prf, max_range)
    
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
  range_bins = end_index - start_index + 1
  
  print("Configuration parameters loaded.")

  return time_stamp, switch_mode, n_seconds, presummed_prf, sampling_rate, range_bins, max_range


def load_data(root_directory, time_stamp, switch_mode, n_seconds, presummed_prf, sample_range):
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
  n_samples_per_channel = n_seconds*presummed_prf*sample_range

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
    n_seconds*presummed_prf, sample_range))
  channel_b_data = np.transpose(channel_b_data.reshape(
    n_seconds*presummed_prf, sample_range))
  
  # Pack data into a matrix
  data = []
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


def rti(data, root_directory, time_stamp, switch_mode, n_seconds, sample_range, presummed_prf, max_range):
  '''
  Produces and saves a range time intensity plot of the data
  Compute a 2D Fast Fourier Transform to convert fast-time axis to range and slow-time axis to time
  '''
  #TODO: Reduce dynamic range of image in order to see more low intensity clutter returns

  # Compute the axes
  x_axis = np.linspace(0, n_seconds, (int)(n_seconds*presummed_prf) , endpoint=False)  # time axis
  y_axis = np.linspace(0, max_range, sample_range, endpoint=False)  # range axis

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
  
  plt.imshow(data_log10, aspect='auto', origin='lower', extent=[min(x_axis), max(x_axis), min(y_axis), max(y_axis)])
    
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


if __name__ == "__main__":
  main()  # call the main function