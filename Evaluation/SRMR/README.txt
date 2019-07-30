README and USER AGREEMENT for the SRMR Matlab Toolbox

This Matlab toolbox computes the adaptive speech-to-reverberation modulation energy
ratio measure (SRMR), as described in the reference below. This version of the program 
operates on 16kHz sampled *speech* files (algorithm has not been tested for audio/music
files). Files using different sampling rates will be automatically converted to 16 kHz.
A simple energy thresholding VAD algorithm is used to remove silence segments longer
than 50ms.

The SRMR toolbox is being provided "as is" under the condition of sole scientific,
non-commercial use. By using the toolbox, you agree to the following definitions and
conditions:

Definitions:
1. User: The person or organization that downloads the SRMR toolbox or any part of it.
2. Provider: Tiago H. Falk, INRS-EMT, Montreal, QC, Canada.

Conditions:
1. The SRMR toolbox is provided without any guarantee.
2. No legal claims of any kind can be made from accepting or using the toolbox.
3. The toolbox provider is not liable for any damage that may result from downloading,
installing, or running the SRMR toolbox.
4. Use of the toolbox is solely for scientific, non-commercial purpose.
5. The provider retains all rights, including copyright and intellectual property
ownership, embodied in the SRMR toolbox.
6. The user will inform the provider of any bugs/errors encountered while using the
toolbox.
7. For any publications which report results obtained using the toolbox, the following
citation should be used:
T. Falk, C. Zheng, W.-Y. Chan, A Non-Intrusive Quality and Intelligibility Measure of
Reverberant and Dereverberated Speech, IEEE Trans Audio Speech Lang Process, Vol.18,
No.7, pp.1766-1774, Sept.2010.

Acknowledgements/Requirements:
Thanks to Mrs. Lu Huo for providing the ITU-T P.56 speech voltmeter. The gammatone
filterbank design uses Malcolm Slaney's Auditory toolbox (the necessary functions have
been embedded).

This script is compatible with MATLAB R2009a or newer. Note that the Matlab Signal Processing 
Toolbox may be needed in order to run the main SRMR script. The script has not been optimized
for processing speed.

Please send comments to tiago.falk@ieee.org.

This version was developed and packaged by João Felipe Santos (jfsantos@emt.inrs.ca)
at June 26th 2013.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Usage: SRMRval=SRMR_main(wavefile);

Sanity check:
The directory "audio" contains a clean speech file F1_010.wav and its reverberant
counterparts with reverberation times (T60) ranging from 0.4-1s (increments of 0.1s)
and 1.5 and 2s. The script "SRMR_test.m" computes SRMR* for all files described
in "test_filenames.txt". The computed SRMR* values for these 10 test files are shown
in "Test_results.xls". The reverberation-to-speech modulation energy ratio (RSMR,
i.e., the inverse of SRMR) is shown to attain a correlation of 0.994 with the true
T60.

