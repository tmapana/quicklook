import configparser
import pandas as pd
import numpy as np
import os.path
import math
import glob
import json
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
        motion_df['abs_vel'].values)
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

def get_rf_freq(FN):
    '''
    Returns the RF frequency associated with a given fractional numerator
    '''
    return RP_CLK*(N_COUNTER + FN/FD)/RF_DIV