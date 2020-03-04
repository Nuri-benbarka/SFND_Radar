# SFND_Radar

The code is organized as follows:
1. FMCW Waveform Design: The FMCW waveform specifications; the chirp time, the
bandwidth and the slope were calculated acording to system specifications which are:
a. The Frequency of operation = 77GHz
b. Max Range = 200m
c. Range Resolution = 1 m
d. Max Velocity = 100 m/s
2. FMCW Signal Generation: The distance of the target was calculated for each time step,
and using the distance the time of travel was calculated. The transmitted signal was
simulated using the mathimatical function of a chrimp function and recieved signal using
its delayed version. Mixing was done by simple element-wise multiplication.
3. Range FFT: An FFT is done on the mixed signal and plotted. The x-axis of the plot is the
range which was calculated using the range estimation formula.
4. 2D CFAR: the CFAR was implemended using the following steps:
a. The test cell is slided across the whole range-doppler matrix but the edges.
b. All the values in the training cells are summed up after converting them into
watts.
c. The threshold is calculated by adding the average of the training cell values and
an offset after converting into dBs again.
d. The test cell is compared to the threshold if it is less than the threshold it is set to
zero otherwise it is set to one.
Selection of Training, Guard cells and offset:
1. Training cells: if a low number of training cells was selected, the threshold will be
less representative of the signal noise. And if a high number of training cells was
selected, the algorithim’s behavior will resemble more the constant threshold
algorithim which will lead to false positives in high noise areas.
2. Guard cells: if no Guard cells are used the signal itself can increase the threshold
then it will cancel itself. And using high number of graud cell will make the
threshold less representative of the signal noise.
3. Offset: if the offset is high everything will be filtered and if set low nothing will be
filtered.
