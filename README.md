# Quicklook SAR Processor

### Base operation:
Read raw data samples and produce IQ matrices for each polarization.
#
### Usage: processor.py dataset_root_directory [OPTIONS]
        OPTIONS:
        'raw' - Save .raw file for observing data more accurately
        'rd' - Plot range-Doppler map for each polarization
        'rti' - Plot and save range time intensity plot for each polarization
        'hist' - Plot histogram of data
        'motion' - Perform motion compensation on the data