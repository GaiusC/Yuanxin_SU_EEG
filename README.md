# Removal of Ocular Artifacts from EEG

## Environments Configuration

### Install TFA toolboxs
Run \Environment\TFA_Toolboxs\install_all_toolboxs.m

### Add toolboxs of EMD and EEMD and CEEMDAN and add the toolbox of FastICA
add path Environments folder and subfolders

## Main
### Descripton
#### LowPass_4Hz_Hamming.m
Function of Lowpass Filter with cur-off frequency 4Hz and Hamming window
#### PSD_Thesholding.m
Function of calculating the powerspectral entropy and information entropy and thresholding processing
#### ICA_OD.m
ICA algorithm and Outlier Detection

## Start
Run \main\Filter_EMD_OD.m
