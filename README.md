# Quicklook SAR Processor

### Base operation:
Read raw data samples and produce IQ matrices for each polarization.
#
### Usage: processor.py dataset_root_directory [OPTIONS]
        OPTIONS:
        'rti' - Plot and save range time intensity plot for each polarization
        'mocomp' - Perform motion compensation on the data
        'sar' - Run G2 Processor to produce SAR image