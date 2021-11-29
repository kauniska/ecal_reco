# ecal_reco
Repository for the ECAL data readout, reconstruction and analysis


I)  Acquirung data with tofpet readout
  1. Start the board-server on DAQ board: 
     * Open a new terminal window
     * Establish a ssh conection with a board
     * Run ```board-server```
  2. Calibrate TOFPET ASICs (if it was done already, skip this step) :
     * In a new terminal window navigate to ```/home/cholak/Software/daq-test```
     * Set HV of SiPMs to 80% of BD voltage
     * Run ```python3 calibration_single.py```
     * Set HV of SiPMs to OP voltage when it's required in the prompt
     * Examine DCR scans using ```plot_dcr.py``` after the calibration process is done
  3. Start the daq-server on host machine: 
     * In a new terminal window navigate to ```/home/cholak/Software/snd_scifi/build-host/bin```
     * Run ```./daq-server ```
  4. Acquire data
     * In a new terminal window navigate to ```/home/cholak/Software/daq-test```
     * Prepare setup, check light-tightness and set HV of SiPMs to OP voltage
     * Run data acquisition using ```test_daq.py``` script

II) Data convertation using calibration files
Acuired data is in the raw format and needs to be converted in order to opbtain physical values. More info can be found [here](https://gitlab.cern.ch/snd-scifi/software/-/wikis/Raw-data-format) and at the [TOFPET data-sheet](/home/cholak/Documents/PETsys TOFPET 2C ASIC - Datasheet (rev 10).pdf)

  1. Raw data is located in ```/home/cholak/Documents/data_test``` and calibrated in ```/home/cholak/Documents/data_conv```
  2. Source data converation utilities with ```source /home/cholak/Documents/data_conv/setup.sh```
  3. Calibration functions can be expoted from: 
