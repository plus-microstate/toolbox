# +microstate toolbox


## About +microstate

+microstate is a MATLAB toolbox for performing microstate analysis in sensor- or source-reconstructed EEG or MEG data. +microstate can be used for individual- or group-level analyses, and applied to both resting-state or task-evoked data. As well as microstate analysis, +microstate includes functions for calculating microstate statistics and microstate-segmented functional connectivity, visualisation, and simulation of data. Basic preprocessing is also available in the toolbox, and for more advanced preprocessing +microstate integrates with a number of pre-exisiting toolboxes for M/EEG processing and analysis. 
## Download and Installation

The easiest way to download the toolbox is to click the following link: 

### [**Download +microstate**](https://github.com/plus-microstate/toolbox/archive/refs/heads/master.zip)

Alternatively, the toolbox can be downloaded from [the toolbox GitHub repository](https://github.com/plus-microstate/toolbox). MATLAB must be installed to run +microstate, which has been tested on MATLAB R2017b and higher. Additionally, the following MATLAB toolboxes are required: [Statistics and Machine Learning](https://www.mathworks.com/products/statistics.html), [Signal Processing](https://www.mathworks.com/products/signal.html), and [Wavelet](https://www.mathworks.com/products/wavelet.html). 

To install the toolbox, unzip the downloaded folder (toolbox-master), then open MATLAB and enter 
```markdown
addpath('PATH/toolbox-master')
```
in the MATLAB console (here `PATH` is the location in which you saved the toolbox-master directory, usually `~/Downloads`). You must either do this each time you open MATLAB, or to permenantly add +microstate to the path you can add toolbox-master to the path using MATLAB's `pathtool` function. 

You can check that the toolbox is properly installed by typing 
```markdown
microstate.functions.toolbox_path
```
in the MATLAB console. If installed correctly, this command returns a string containing the directory with the +microstate toolbox. Otherwise, you will receive an error.

The first time you run the toolbox, you need to install some additional files. This is done by typing
```
microstate.functions.install
```
This does not need to be run every time you use +microstate, but it is worth running occasionally since the command will additionally update the toolbox to the latest version if an older version is installed. 

## Getting started

For your first microstate analysis using +microstate, we recommend working through the tutorials. The toolbox-master folder contains a subfolder called "examples", which contains MATLAB Live Scripts to perform the four example analyses described in the [toolbox manuscript](https://doi.org/10.1101/2021.07.13.452193). We recommend working through each tutorial in order, as these analyses increase in complexity and assume knowledge from previous examples. 

For extra help, you can refer to the following resources: 
- [Toolbox manuscript](https://doi.org/10.1101/2021.07.13.452193)
- [Wiki pages](https://github.com/plus-microstate/toolbox/wiki)
- [Our original paper](https://doi.org/10.1101/2021.03.25.436979) - particularly Supplementary Material - for methodological details 

## Citing the toolbox
If you use the +microstate toolbox in your analysis, please cite the toolbox paper (reference given below) and include the URL to the toolbox webpage. 

**Example citation**: 

Microstate segmentation, analysis, and visualisation used the +microstate toolbox for brain microstate analysis in sensor and cortical EEG/MEG (Tait and Zhang (2021); [plus-microstate.github.io](http://plus-microstate.github.io)). 

**Reference**:

Tait and Zhang (2021), _+microstate: A MATLAB toolbox for brain microstate analysis in sensor and cortical EEG/MEG_, bioRxiv 2021.07.13.452193. doi: [10.1101/2021.07.13.452193](https://doi.org/10.1101/2021.07.13.452193)

## Bug reports and contact

If you encounter and issues or bugs, please contact us via one of the following means: 
- [GitHub Issue Tracker](https://github.com/plus-microstate/toolbox/issues)
- [Email](mailto:l.tait@bham.ac.uk)

For queries about use of the toolbox, feedback, suggestions, or interest in collaboration, please contact [Luke Tait via email](mailto:l.tait@bham.ac.uk). There is additionally an informal mailing list, where we will occasionally send out information about toolbox updates or relevant publications. If you are interested in joining this mailing list please use the link above to contact Luke Tait via email. 

## Publications

- Tait and Zhang (2022), _MEG cortical microstates: spatiotemporal characteristics, dynamic functional connectivity and stimulus-evoked responses_, NeuroImage 251:119006. doi: [10.1016/j.neuroimage.2022.119006](https://doi.org/10.1016/j.neuroimage.2022.119006). 
Paper describing the methods underpinning the toolbox. 
- Tait and Zhang (2021), _+microstate: A MATLAB toolbox for brain microstate analysis in sensor and cortical EEG/MEG_, bioRxiv 2021.07.13.452193. doi: [10.1101/2021.07.13.452193](https://doi.org/10.1101/2021.07.13.452193). 
Preprint describing the format of the toolbox and tutorials. 
- Tait et al. (2020) _EEG microstate complexity for aiding early diagnosis of Alzheimer’s disease_, Sci Rep 10:17627. doi: [10.1038/s41598-020-74790-7](https://doi.org/10.1038/s41598-020-74790-7). 
Early publication using many of the codes which later became the +microstate toolbox to perform EEG microstate analysis in patients with Alzheimer's disease. 


