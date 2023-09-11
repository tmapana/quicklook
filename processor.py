"""
This is a quicklook processor which is used to give a quick
snapshot of miloSAR data for review

Author:   Tlotliso Mapana
Date:     28 September 2021
"""

# TODO: Motion compensation + clean up G2

import os
import sys
import glob
import json
import numpy as np
import configparser
import matplotlib.pyplot as plt

from string import *
from PIL import Image
from matplotlib import cm
from scipy import interpolate

import g2tools
from parameters import *

# TODO: Motion compensation
# TODO: Basic SAR processor
# TODO: Feed-through filter

#---------------------------------------------------------------------------------------------------------------------------------#
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
    time_stamp, switch_mode, n_seconds, prf, decimation_factor, presummed_prf, ns_pri, \
        min_range, max_range, n_range_bins, sampling_rate, g2_adc, range_resolution, \
          ns_fft, centre_freq, wavelength, n_chunks \
              = load_radar_parameters(root_directory=root_directory)
    
    # Axes definitions
    ramp_period = 1/prf
    chirp_rate = ramp_period/RAMP_BANDWIDTH
    range_scaling = (c/2) / chirp_rate

    # range_axis = np.linspace(0, sampling_rate, ns_fft, endpoint=False) * range_scaling
    # time_axis = np.linspace(0, n_seconds, int(n_chunks), endpoint=False)
    frequency_axis = np.linspace(-sampling_rate/2, sampling_rate/2, ns_fft, endpoint=False)

    # load data from file and pack into a data matrix
    # also plots signal envelope for analysis
    # option to remove RFI by replacing range bins containing RFI with zero arrays
    data = load_data(root_directory=root_directory, time_stamp=time_stamp, switch_mode=switch_mode, \
                    n_seconds=n_seconds, presummed_prf=presummed_prf, ns_pri=ns_pri, check_rfi=True)

    data = fast_time_fft(data=data, ns_fft=ns_fft)

    power_spectrum(data=data, root_directory=root_directory, ns_fft=ns_fft, frequency_axis=frequency_axis, n_chunks=n_chunks)

    # OPTIONS
    if len(sys.argv) > 2:
      for opt in range(2, len(sys.argv)):
        if sys.argv[opt] == 'raw':
          # save binary image to file
          save_raw(data, root_directory, switch_mode, n_range_bins, n_seconds, time_stamp, ns_fft)

        elif sys.argv[opt] == 'notch':
          # data = notch_filter(data=data, sampling_rate=5e9, intereference_freq=2.4e9, notch_width=20e6)
          data = notch_filter(data=data, n_chunks=n_chunks, intergerence_frequency=2.4e9, quality_factor=30, sampling_rate=sampling_rate)

        elif sys.argv[opt] == 'cfar':
          # execute constant false alarm rate filter tp supress RFI
          # data = cfar(root_directory=root_directory, data=data, ns_fft=ns_fft, switch_mode=switch_mode, \
                # n_chunks=n_chunks, ns_cpi=n_chunks)
          
          data = cfar_detector(root_directory=root_directory, data=data, ns_fft=ns_fft, n_chunks=n_chunks, \
            n_guard=2, n_train=10, pfa=1e-3)

        elif sys.argv[opt] == 'rd':
          # produce a range-Doppler plot and save to directory
          range_doppler(data, root_directory, time_stamp, switch_mode, \
                        n_seconds, presummed_prf, max_range, sampling_rate, all=False)

        elif sys.argv[opt] == 'rti':
          # produce range-time plot and save to diectory
          rti(data, root_directory, time_stamp, switch_mode, \
              n_seconds, ns_pri, presummed_prf, max_range, title='rti', all=False)

        elif sys.argv[opt] == 'rtp':
          # produce range-time plot and save to diectory
          rti(data, root_directory, time_stamp, switch_mode, \
              n_seconds, ns_pri, presummed_prf, max_range, title='rtp', all=False)

        elif sys.argv[opt] == 'mocomp':
          # perform motion compensation
          data = motion_compensation(data=data, switch_mode=switch_mode, root_directory=root_directory, \
                                      n_seconds=n_seconds, n_range_bins=int(ns_fft), n_chunks=n_chunks, \
                                        wavelength=wavelength, max_range=max_range, prf=prf, \
                                          range_resolution=range_resolution, ns_fft=ns_fft, plot_path=True)
        
        elif sys.argv[opt] == 'sar':
          # run SAR processor
          SAR(data=data, root_directory=root_directory, switch_mode=switch_mode, time_stamp=time_stamp, \
              prf=presummed_prf, n_range_bins=n_range_bins, max_range=max_range, g2_adc=g2_adc, \
              range_resolution=range_resolution, n_chunks=n_chunks, centre_freq=centre_freq, d_range=[])

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

#---------------------------------------------------------------------------------------------------------------------------------#
def load_data(root_directory, time_stamp, switch_mode, n_seconds, presummed_prf, ns_pri, check_rfi=False):
  '''
  This function reads in raw miloSAR data stored in a .bin file in the given directory
  '''
  # open file stored at existing directory and read as binary
  file_name = glob.glob(os.path.join(root_directory, '*.bin'))[0]  
  with open(file_name, 'rb') as f:
    # read data into an array as signed 16-bit integer values
    data_array = np.fromfile(f, np.int16)
  
  print("\nLoading dataset...")
  
  '''
  data_array is an array of interleaved IQ data samples.
  Each sample contains I and Q data point from two channels, A and B
  Need to extract each channel data respectively and reshape data into a 2D matrix
  '''
  # each data sample contains [Q_a, I_a, Q_b, I_b]
  channel_a_data = data_array[1::4] + 1j*data_array[0::4]
  channel_b_data = data_array[3::4] + 1j*data_array[2::4]

  # calculate the number of samples recorded per channel
  n_chunks = int(n_seconds*presummed_prf)
  n_samples_per_channel = int(n_chunks*ns_pri)

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
    n_chunks, ns_pri))
  channel_b_data = np.transpose(channel_b_data.reshape(
    n_chunks, ns_pri))
  
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

  # scale data
  data = data*ADC_BITS_VOLTS

  # apply windowing function to data
  data = window(data=data, switch_mode=switch_mode, ns_pri=ns_pri)

  # attempt to remove range bins contamited by RFI
  if check_rfi:
    remove_rfi(data=data, ns_pri=ns_pri, n_chunks=n_chunks)

  # view signal amplitude in the time-domain
  signal_envelope = pow(abs(data[:, 0:1000, 1]), 2)
  signal_envelope = np.sum(signal_envelope, axis=1)
  signal_envelope = 2*linear2db(signal_envelope)

  plt.figure()
  plt.plot(signal_envelope)
  plt.xlim(0, ns_pri)
  plt.ylim(-100, 0)
  plt.grid()
  plt.xlabel('Range Bins')
  plt.ylabel('Amplitude [dB]')
  plt.title('Signal Envelope (1000 Pulses)')
  plt.legend(['ChA', 'ChB'])
  if check_rfi:
    plt.savefig(os.path.join(root_directory, 'quicklook/signal_envelope_rfi_removed.svg'))
  else:
    plt.savefig(os.path.join(root_directory, 'quicklook/signal_envelope.svg'))
  
  # data_cfar = cfar_filter(data=data, n_training_cells=10000)

  print("Dataset", time_stamp, "loaded successfully.")
  return data

#---------------------------------------------------------------------------------------------------------------------------------#
def load_radar_parameters(root_directory):
    '''
    Important information about the data is stored in the configuration file of each data directory
    This function reads in the config file and extracts the following experiment parameters:
    1. switch_mode = integer indicating which channel(s) was used during operation
    2. n_seconds = length of data take in seconds
    3. prf = pulse repetition frequency
    4. sampling_rate = sampling rate of ADC
    5. ns_pris = factor applied to the prf when writing data to physical disk
    6. start_index = position to start recording data
    7. end_index = position to end recording
    '''
    summary = configparser.ConfigParser()
    summary.read(os.path.join(root_directory, "summary.ini"))

    time_stamp = summary['general']['time_stamp']
    switch_mode = int(summary['dataset']['switch_mode'])
    n_seconds = int(summary['dataset']['n_seconds'])
    prf = int(summary['dataset']['prf'])
    decimation_factor = int(summary['dataset']['decimation_factor'])
    sampling_rate = float(summary['dataset']['sampling_rate'])
    presumming_factor = int(summary['integration']['n_pris'])
    start_index = int(summary['integration']['start_index'])
    end_index = int(summary['integration']['end_index'])

    setup = configparser.ConfigParser()
    setup.read(os.path.join(root_directory, "setup.ini"))

    min_range_string = setup['geometry']['min_range']  # 100.0; [m]
    min_range = int(min_range_string.split('.')[0])
    max_range_string = setup['geometry']['max_range']  # 900.0; [m]
    max_range = int(max_range_string.split('.')[0])

    print("Configuration parameters loaded.")   # check point

    # realisable prf of data file stored in memory as data rate is halved when writing to memory
    presummed_prf = prf/presumming_factor # TODO rename to more understandable term

    # number of range bins recorded for each pri
    ns_pri = end_index - start_index + 1

    tx_fn_init = int(summary['tx_synth']['fractional_numerator'])
    ramp_increment = int(summary['tx_synth']['up_ramp_increment'])
    ramp_length = int(summary['tx_synth']['up_ramp_length'])
    ramp_period = ramp_length/RP_CLK

    ramp_start_freq = get_rf_freq(tx_fn_init)
    ramp_end_freq = get_rf_freq(ramp_increment*ramp_length + tx_fn_init)
    ramp_bandwidth = ramp_end_freq - ramp_start_freq

    chirp_rate = ramp_bandwidth/ramp_period

    centre_freq = (ramp_start_freq + ramp_end_freq)/2
    range_scaling = (c/2)/chirp_rate
    wavelength = c/centre_freq

    ns_fft = next_pow_two(ns_pri)

    spectral_resolution = sampling_rate/ns_pri # 1/T
    fft_bin_spacing = sampling_rate/ns_fft
    range_bin_spacing = fft_bin_spacing * range_scaling
    range_resolution = range_scaling * np.hypot(spectral_resolution, fft_bin_spacing)
    g2_adc = c / (2*range_bin_spacing)

    min_range_bin = int(np.floor(min_range / range_bin_spacing))
    max_range_bin = int(np.ceil(max_range / range_bin_spacing))

    n_range_bins = max_range_bin - min_range_bin

    min_chunk = 0
    max_chunk = int(np.ceil(presummed_prf*n_seconds))
    n_chunks = max_chunk - min_chunk
    if switch_mode == 3:
        n_chunks = int(n_seconds*presummed_prf//2)

    return time_stamp, switch_mode, n_seconds, prf, decimation_factor, presummed_prf, ns_pri, \
            min_range, max_range, n_range_bins, sampling_rate, g2_adc, range_resolution, \
              ns_fft, centre_freq, wavelength, n_chunks

#TODO ---------------------------------------------------------------------------------------------------------------------------------#
def save_raw(data, root_directory, switch_mode, n_range_bins, n_seconds, time_stamp):
  '''
  Not currently saving raw data as intended
  TODO: Rewrite this function to produce correct results
  '''
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

#TODO ---------------------------------------------------------------------------------------------------------------------------------#
def range_doppler(data, root_directory, time_stamp, switch_mode, n_seconds, presummed_prf, max_range, sampling_rate, all=False):
  '''
  Produces and saves a Range_Doppler plot of SAR data
  Convert fast-time (axis=1) to range and slow-time (axis=0) to Doppler by way
  of the Fast Fourier Transform
  '''

  # Set FFT length
  n_FFT_Doppler = 2048

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

#---------------------------------------------------------------------------------------------------------------------------------#
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
  data_min = data_max - 15
  data_dB = np.clip(data_dB, data_min, data_max)
    
  plt.imshow(X=data_dB, aspect='auto', origin='lower', extent=[min(x_axis), max(x_axis), min(y_axis), max(y_axis)])
  plt.colorbar(label='Power [dB]')
  plt.tight_layout()
    
  file_name = time_stamp + "_" + title + "_" + str(switch_mode) + channel + ".png"
  plt.savefig(os.path.join(root_directory, 'quicklook/'+file_name), dpi=300)
  
  print("\n" + time_stamp + " " + title.upper() + " " + str(switch_mode) + channel + " saved.")

#---------------------------------------------------------------------------------------------------------------------------------#
def rti(data, root_directory, time_stamp, switch_mode, n_seconds, ns_pri, presummed_prf, max_range, title, all=False):
  '''
  Produces and saves a range time intensity plot of the data
  Compute a 2D Fast Fourier Transform to convert fast-time axis to range and slow-time axis to time
  '''
  #TODO: Reduce dynamic range of image in order to see more low intensity clutter returns

  # Compute the axes
  x_axis = np.linspace(0, n_seconds, ns_pri, endpoint=False)  # time axis
  y_axis = np.linspace(0, max_range, (int)(n_seconds*presummed_prf), endpoint=False)  # range axis

  # select which channel to plot depending on switch mode
  if switch_mode == 1:
    # unpack data into the relevent polarizations depending on hardware configuration
    channel_1a = data[:,:,0]
    channel_1b = data[:,:,1]

    # channel_1a_fft = fast_time_fft(data=channel_1a, ns_fft=ns_fft)
    plot_rti(title=title, data=channel_1a, x_axis=x_axis, y_axis=y_axis, switch_mode=switch_mode, \
              channel='a', root_directory=root_directory, time_stamp=time_stamp)

    if all:
      plot_rti(title=title, data=channel_1b, x_axis=x_axis, y_axis=y_axis, switch_mode=switch_mode, \
              channel='b', root_directory=root_directory, time_stamp=time_stamp)
  
  elif switch_mode == 2:
    # unpack data into the relevent polarizations depending on hardware configuration
    channel_2a = data[:,:,0]
    channel_2b = data[:,:,1]

    # FFT in the fast-time axis
    plot_rti(title=title, data=channel_2a, x_axis=x_axis, y_axis=y_axis, switch_mode=switch_mode, \
              channel='a', root_directory=root_directory, time_stamp=time_stamp)

    if all:
      # FFT in the fast-time axis
      plot_rti(title=title, data=channel_2b, x_axis=x_axis, y_axis=y_axis, switch_mode=switch_mode, \
              channel='b', root_directory=root_directory, time_stamp=time_stamp)

  elif switch_mode == 3:
    channel_1a = data[:,:,0]
    channel_1b = data[:,:,1]
    channel_2a = data[:,:,2]
    channel_2b = data[:,:,3]
    
    # channel_1a_fft = fast_time_fft(data=channel_1a, ns_fft=ns_fft)
    plot_rti(title=title, data=channel_1a, x_axis=x_axis, y_axis=y_axis, switch_mode=1, \
              channel='a', root_directory=root_directory, time_stamp=time_stamp)
    
    if all:      
      plot_rti(title=title, data=channel_1b, x_axis=x_axis, y_axis=y_axis, switch_mode=1, \
              channel='b', root_directory=root_directory, time_stamp=time_stamp)
          
      plot_rti(title=title, data=channel_2a, x_axis=x_axis, y_axis=y_axis, switch_mode=2, \
              channel='a', root_directory=root_directory, time_stamp=time_stamp)
      
      plot_rti(title=title, data=channel_2b, x_axis=x_axis, y_axis=y_axis, switch_mode=2, \
              channel='b', root_directory=root_directory, time_stamp=time_stamp)

#---------------------------------------------------------------------------------------------------------------------------------#
def plot_rti(title, data, x_axis, y_axis, switch_mode, channel, root_directory, time_stamp):
  '''
  Plot and save Range-Time Intensity to file
  '''
  plt.subplots(nrows=1,ncols=1)
  plt.title(title.upper() + " " + str(switch_mode) + channel)
  plt.xlabel("Slow Time [s]")
  plt.ylabel("Range [m]")
  
  if title == 'rti':
    data = data/np.amax(data)
    data_dB = 20*np.log10(abs(data))
    data_max = np.amax(data_dB)
    data_min = np.nanmin(np.mean(data_dB, axis=1)) 

    data_dB = np.clip(data_dB, data_min, data_max)

    plt.imshow(
      data_dB, 
      interpolation='none', 
      cmap='viridis', 
      aspect='auto', 
      origin='lower', 
      extent=[min(x_axis), max(x_axis), min(y_axis), max(y_axis)])
  
  elif title == 'rtp':
    data_max = np.pi
    data_min = -np.pi

    plt.imshow(
      np.angle(data),
      interpolation='none',
      aspect='auto', 
      origin='lower', 
      vmax=data_max,
      vmin=data_min,
      extent=[min(x_axis), max(x_axis), min(y_axis), max(y_axis)])
      
  file_name = time_stamp + "_" + title + "_" + str(switch_mode) + channel + ".png"
  plt.savefig(os.path.join(root_directory, 'quicklook/'+file_name), dpi=300)
  
  print("\n" + time_stamp + " " + title.upper() + " " + str(switch_mode) + channel + " saved.")

#---------------------------------------------------------------------------------------------------------------------------------#
def motion_compensation(data, switch_mode, root_directory, n_seconds, n_range_bins, n_chunks, wavelength, \
                        max_range, prf, range_resolution, ns_fft, plot_path=True):
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
  print("\nCorrecting motion errors...")
  x_axis = np.linspace(0, n_seconds, n_range_bins, endpoint=False)  # time axis
  y_axis = np.linspace(0, max_range, n_chunks, endpoint=False)  # range axis
  motion_filename = glob.glob(os.path.join(root_directory, '*.json'))[0]
  
  # read motion file contents
  motion_data = [json.loads(line) for line in open(motion_filename, 'r')]

  timestamp = np.array([])
  latitude = np.array([])
  longitude = np.array([])
  altitude = np.array([])

  for line in motion_data:
    timestamp = np.append(timestamp, line['timestamp'])
    latitude = np.append(latitude, line['lat'])
    longitude = np.append(longitude, line['lon'])
    altitude = np.append(altitude, line['u'])
  
  # form nominal flight path from start and end lat-long coordinates - straight line fit
  fit = np.polyfit(np.array([longitude[0], longitude[-1]]), np.array([latitude[0], latitude[-1]]), deg=1)
  
  # nominal path will be a line through (longitude, nom_path(longitude))
  nom_path = np.poly1d(fit) # the nominal path the plane should have taken
  
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

  range_axis = np.linspace(0, n_seconds, n_range_bins, endpoint=False)
  range_deviation = interpolate.interp1d(t_interpolate, range_deviation, kind='linear')
  range_deviation = range_deviation(range_axis)
  
  # apply phase correction
  # TODO: extract system hardware properties correctly
  phase_shift = np.exp(np.multiply(-1j*(4*np.pi)/wavelength, range_deviation))
  # phase_shift = np.dstack((phase_shift,))

  data_out = data
  # data_out = np.multiply(data, np.transpose(phase_shift))
  for rng in range(0,n_range_bins):
    data_out[rng,:,:] = np.multiply(data[rng,:,:], phase_shift[rng])

  # apply range bin correction
  #   find max r_fft
  #   determine range bin size = r_fft/ns_fft
  # total_range = max_range - min_range
  range_bin_size = max_range/n_range_bins # [m] each range bin spans

  # number of range bins to shift given delta_range calculated above
  range_bin_shift = np.round(range_deviation/range_bin_size)
  
  if switch_mode==1 or switch_mode==2:
    temp = np.zeros((n_range_bins, n_chunks, 2)).astype('complex64')
  elif switch_mode==3:
    temp = np.zeros((n_range_bins, n_chunks, 4)).astype('complex64')

  count_temp = 0
  # for i in range(0, n_chunks):
  for j in range(0, n_range_bins):
    x = int(range_bin_shift[j])
    
    # limit cases where range shift is too large
    # if x > 3:
      # x = 3
    # elif x < -3:
      # x = -3
    
    new_rbin = j - int(x)

    # shift range bin from j to new_rbin
    if (new_rbin >= 0 and new_rbin < n_range_bins):
      count_temp += 1
      temp[j, :] = data_out[new_rbin, :]

  # data_out = temp
  print(count_temp)

  '''print("\nExtracting target height information form scene...")

  altitude_interp = interpolate.interp1d(t_interpolate, altitude, kind='linear')
  altitude_interp = altitude_interp(range_axis)
  max_alt = min(altitude_interp)

  target_height = np.zeros((n_range_bins, n_chunks))
  elevation_angle = np.deg2rad(90-DEPRESSION_ANGLE)
  data_peak = data/np.amax(np.abs(data))
  for t in range(0, n_range_bins):
    t_height = max_alt - np.cos(elevation_angle)*20*np.log10(abs(data[t,:,0]))*range_resolution
    # slant_range = abs(data[t,:,0]*(c/(2*prf)))
    # t_height = altitude_interp[t] - np.cos(elevation_angle)*slant_range
    target_height[t] = t_height
  
  fig = plt.figure()
  ax = plt.axes(projection='3d')
  X, Y = np.meshgrid(x_axis, y_axis)

  # platform path
  ax.plot_surface(X, Y, np.transpose(target_height), linewidth=0, antialiased=False)
  fig.savefig(os.path.join(root_directory, 'quicklook/scene_height.png'))

  # print("altitude:",altitude_interp[100])
  # print("taget height",target_height[100][1000])'''

  print("Motion errors corrected on dataset.")
  return temp

#---------------------------------------------------------------------------------------------------------------------------------#
def SAR(data, root_directory, switch_mode, time_stamp, prf, n_range_bins, max_range, g2_adc, range_resolution, \
        n_chunks, centre_freq, d_range=[]):
  '''
  Make call to G2 program
  '''
  print("\nRunning the SAR processor...")

  # calculate synthetic aperture from motion data
  synthetic_aperture = load_motion_data(root_directory=root_directory, parameter='aperture')

  # pixel scaling for SAR image
  scaled_width = int(n_range_bins / max_range * synthetic_aperture)
  scaled_height = n_range_bins
  G2_out = data

  if switch_mode == 1:
    G2_out_1a = G2(data=data[:,:,0], title='1a', root_directory=root_directory, time_stamp=time_stamp, \
            prf=prf, n_range_bins=n_range_bins, g2_adc=g2_adc, range_resolution=range_resolution, n_chunks=n_chunks, centre_freq=centre_freq)
    G2_out_1b = G2(data=data[:,:,1], title='1b', root_directory=root_directory, time_stamp=time_stamp, \
            prf=prf, n_range_bins=n_range_bins, g2_adc=g2_adc, range_resolution=range_resolution, n_chunks=n_chunks, centre_freq=centre_freq)
    G2_out = np.dstack((G2_out_1a, G2_out_1b))

  elif switch_mode == 2:
    G2_out_2a = G2(data=data[:,:,0], title='2a', root_directory=root_directory, time_stamp=time_stamp, \
            prf=prf, n_range_bins=n_range_bins, g2_adc=g2_adc, range_resolution=range_resolution, n_chunks=n_chunks, centre_freq=centre_freq)
    G2_out_2b = G2(data=data[:,:,1], title='2b', root_directory=root_directory, time_stamp=time_stamp, \
            prf=prf, n_range_bins=n_range_bins, g2_adc=g2_adc, range_resolution=range_resolution, n_chunks=n_chunks, centre_freq=centre_freq)

    G2_out = np.dstack((G2_out_2a, G2_out_2b))
  
  elif switch_mode == 3:
    
    G2_out_1a = G2(data=data[:,:,0], title='1a', root_directory=root_directory, time_stamp=time_stamp, \
            prf=prf, n_range_bins=n_range_bins, g2_adc=g2_adc, range_resolution=range_resolution, n_chunks=n_chunks, centre_freq=centre_freq)
    G2_out_1b = G2(data=data[:,:,1], title='1b', root_directory=root_directory, time_stamp=time_stamp, \
            prf=prf, n_range_bins=n_range_bins, g2_adc=g2_adc, range_resolution=range_resolution, n_chunks=n_chunks, centre_freq=centre_freq)
    G2_out_2a = G2(data=data[:,:,2], title='2a', root_directory=root_directory, time_stamp=time_stamp, \
            prf=prf, n_range_bins=n_range_bins, g2_adc=g2_adc, range_resolution=range_resolution, n_chunks=n_chunks, centre_freq=centre_freq)
    G2_out_2b = G2(data=data[:,:,3], title='2b', root_directory=root_directory, time_stamp=time_stamp, \
            prf=prf, n_range_bins=n_range_bins, g2_adc=g2_adc, range_resolution=range_resolution, n_chunks=n_chunks, centre_freq=centre_freq)
  
    G2_out= np.dstack((G2_out_1a, G2_out_1b))
    G2_out = np.dstack((G2_out, G2_out_2a))
    G2_out = np.dstack((G2_out, G2_out_2b))
  
  # noncoherent avering of data
  image = pow(np.mean(abs(G2_out), axis=2), 2)
  image.astype('complex64').tofile(
			os.path.join(root_directory, 'quicklook/' + 'image.bin') )
  
  # normalize data
  image = image/np.amax(image)
    
  plt.subplots(nrows=1,ncols=1)
  plt.xlabel('Azimuth [m]')
  plt.ylabel('Range [m]')

  image = 20*np.log10(abs(image))
  img_mean = np.mean(image, axis=1)
  
  # dynamic range
  a_min = np.nanmin(img_mean[img_mean!=-np.inf])
  a_max = np.amax(image)
  if d_range:
    if d_range[0]:
      a_min = d_range[0]
    if d_range[1]:
      a_min = d_range[1]
  
  image = np.clip(image, a_min=a_min, a_max=a_max)
  
  # cmap = [viridis, plasma, inferno, magma, cividis, PiYG, PRGn, BrBG, PuOr, RdGy, RdBu, RdYlBu, RdYlGn, Spectral, coolwarm, bwr, seismic]
  plt.imshow(image, cmap='viridis', origin='upper', aspect='equal', interpolation='none', vmin=None, vmax=None)
  file_name = time_stamp + '.sar.noco'
  image_path = os.path.join(root_directory, 'quicklook/'+file_name)
  plt.imsave(image_path, image, cmap='viridis', origin='upper', format='png')#, dpi=300)

  sar_image = Image.open(image_path)
  sar_image.load()
  sar_image = sar_image.resize((scaled_width, scaled_height), resample=Image.LANCZOS)
  sar_image.save(os.path.splitext(image_path)[
               0] + os.path.splitext(image_path)[1], format='png')
    
  print("SAR image saved.")

#---------------------------------------------------------------------------------------------------------------------------------#
def G2(data, title, root_directory, time_stamp, prf, n_range_bins, g2_adc, range_resolution, n_chunks, centre_freq):
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

  mean_velocity = load_motion_data(root_directory=root_directory, parameter='velocity')

  with open(g2_cmd, 'w') as f_id:
    f_id.write('miloSAR azimuth compression command file (azcom)\n')
    f_id.write('$ProgramVersion (jmh)         => 1.1\n\n')

    f_id.write('$ScreenUpdateRate             => ' + str(1) + '\n')
    f_id.write('$LogFileName                  => ' +
            str(g2_log) + '\n')
    f_id.write('$InputStartSampleDelay        => ' + str(0) + '\n')
    f_id.write('$CarrierFreq [Hz]             => ' + str(centre_freq) + '\n')
    f_id.write('$InputPRF [Hz]                => ' + str(prf) + '\n')
    f_id.write('$NomGroundSpeed [m/s]         => ' + str(mean_velocity) + '\n')
    f_id.write('$InputFileAzPts               => ' + str(n_chunks) + '\n')
    f_id.write('$StartProcessAzPt             => ' + str(0) + '\n')
    f_id.write('$AzPtsToProcess               => ' + str(n_chunks) + '\n')
    f_id.write('$InputFileRngBins             => ' + str(n_range_bins) + '\n')
    f_id.write('$StartProcessRngBin           => ' + str(0) + '\n')
    f_id.write('$RngBinsToProcess             => ' + str(n_range_bins) + '\n')
    f_id.write('$InputDCOffsetI               => ' + str(0.0) + '\n')
    f_id.write('$InputDCOffsetQ               => ' + str(0.0) + '\n')
    f_id.write('$InvFFTSizeReduc [pow of 2]   => ' + str(1) + '\n')
    f_id.write('$InputFileName                => ' + str(g2_input) + '\n')
    f_id.write('$OutputFileName               => ' + str(g2_out) + '\n')
    f_id.write('$AppendExistOutFileFlg [Y/N]  => ' + str('N') + '\n')
    f_id.write('$RngFocSegments               => ' + str(-1) + '\n')
    f_id.write('$RefFuncSign [+-1]            => ' + str(1) + '\n') #change
    f_id.write('$A2DFreq [Hz]                 => ' + str(g2_adc) + '\n')
    f_id.write('$NomAzRes [m]                 => ' + str(range_resolution/2) + '\n')
    f_id.write('$WinConstTime [0.0-1.0]       => ' + str(0.08) + '\n')
    f_id.write('$NumLooks                     => ' + str(1) + '\n')
    f_id.write('$LookOverlapFrac [0.0-1.0]    => ' + str(0.0) + '\n')
    f_id.write('$WinConstFreq [0.0-1.0]       => ' + str(0.08) + '\n')
    f_id.write('$RngCurvInterpSize            => ' + str(16) + '\n')
    f_id.write('$RngCurvBatchSize             => ' + str(16) + '\n') #play with this
    f_id.write('$PostSumRatio                 => ' + str(1) + '\n')
    f_id.write('$DetectMethod                 => ' + str(0) + '\n')
    f_id.write('$InputDataType                => ' + str(3) + '\n')
    f_id.write('$OutputDataType               => ' + str(3) + '\n')
    f_id.write('$Scale                        => ' + str(1) + '\n')
    f_id.write('$ReportMax [1/0]              => ' + str(1) + '\n')

  azcom_executable = os.path.join(os.getcwd(), 'g2', 'azcom')
  # os.system(azcom_executable + ' ' + str(g2_cmd) + ' > ' + root_directory + '/quicklook/g2out.txt')
  os.system(azcom_executable + ' ' + str(g2_cmd))

  image = np.fromfile(g2_out, dtype='complex64')
  G2_out = np.flipud(image.reshape(n_range_bins, n_chunks))
  G2_out = np.nan_to_num(G2_out)

  return G2_out

#---------------------------------------------------------------------------------------------------------------------------------#
def range_profile(data, root_directory, range_bin, range_axis):
  '''
  Plot range profile of given range line
  '''
  plt.figure()
  plt.plot(range_axis, linear2db(abs(data[:, range_bin])/np.amax(data)))
  plt.xlim(0, range_axis[-1])
  plt.grid()
  plt.xlabel('Range (m)')
  plt.ylabel('Normalised Power [dB]')
  plt.title('Range Profile of Range Line', range_bin)
  title = 'quicklook/range_profiles/'+range_bin+'.png'
  plt.savefig(os.path.join(root_directory, title))

#---------------------------------------------------------------------------------------------------------------------------------#
def power_spectrum(data, root_directory, ns_fft, frequency_axis, n_chunks):
  '''
  Power spectum of the entire signal
  To determine in which range lines RFI interference is most prevalent
  Store power spectrum values in an array to determine range lines with high returns
  '''  
  power_spectrum = pow(abs(data[:,0:1000,:]), 2)
  power_spectrum = np.sum(power_spectrum, axis=1)/n_chunks
  power_spectrum = linear2db(power_spectrum)

  plt.figure()
  plt.plot(frequency_axis/1e6, power_spectrum)
  # plt.xlim(min(frequency_axis), max(frequency_axis))
  plt.ylim(np.amin(power_spectrum)-10, 0)
  plt.grid()
  plt.xlabel('Frequency [MHz]')
  plt.ylabel('Power [dB]')
  plt.title('Power Spectrum of Signal (1000 Pulses)')
  plt.legend(['ChA', 'ChB'])
  plt.savefig(os.path.join(root_directory, 'quicklook/power_spectrum.svg'))


  # Find the range line with the max 

def remove_rfi(data, ns_pri, n_chunks):

  # find amplitude of each range bin
  # remove range line with presence of rfi

  signal_envelope = pow(abs(data[:, 0:n_chunks,0]), 2)
  signal_envelope = np.average(signal_envelope, axis=1)
  signal_envelope = 2*linear2db(signal_envelope)

  n_channels = data.shape[-1]
  data_clean = data
  # iterate through the data and identify range bins with rfi presence and zero them out
  for rng_bin in range(1, ns_pri):
    prev = signal_envelope[rng_bin-1]
    current = signal_envelope[rng_bin]

    # if (abs(current)-abs(prev) > 3):  # if the current range bin contains twice or more of the power of the previous range bin
      # data_clean[rng_bin,:,:] = np.zeros((n_chunks,n_channels), dtype=np.int16)

  return data_clean
  
#---------------------------------------------------------------------------------------------------------------------------------#
if __name__ == "__main__":
  main()  # call the main function