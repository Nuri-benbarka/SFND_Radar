clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8

range_max = 200;
range_res = 1;
vel_max = 70;
vel_res = 3;
c = 3e8;

%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
R_init = 110;
vel = -20;


%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
Ts = (5.5 * 2 * range_max) / c ;

Bandwidth = c / (2 * range_res);

a = Bandwidth / Ts;


%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq

                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Ts,Nr*Nd); %total time for samples

%distance over time
R = R_init + vel * t;
delay = (2 * R) / c;
Tx = cos(2 * pi * (fc *  t           + ((a *  rem(t,Nr)         .^2) / 2)));
Rx = cos(2 * pi * (fc * (t - delay)  + ((a * (rem(t,Nr) - delay).^2) / 2)));
Mix = Tx .* Rx;


%% RANGE MEASUREMENT


 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix = reshape(Mix,[Nr,Nd]);

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
Mix_fft = fft(Mix,Nr);
Mix_fft = Mix_fft / Nr;

 % *%TODO* :
% Take the absolute value of FFT output
Mix_fft_2s = abs(Mix_fft);

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
Mix_fft_1s = Mix_fft_2s(1:Nr/2 + 1,:);
%Mix_fft_1s(2:end-1,:) = 2 * Mix_fft_1s(2:end-1,:);
f = (1/t(2))*(0:(Nr/2))/Nr;
range = (c * Ts * f) / (2 * Bandwidth);


%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

 % *%TODO* :
 % plot FFT output
plot(range,Mix_fft_1s(:,1));

 
axis ([0 range_max 0 0.5]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map
%Select the number of Training Cells in both the dimensions.
Tr=20;
Td=18;


%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr=3;
Gd=2;


% offset the threshold by SNR value in dB
offset=8;


%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(size(RDM));
Cut = zeros(size(RDM));



%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.
for i=Tr+Gr+1:Nr/2-(Gr+Tr)
    for j = Td+Gd+1:Nd-(Gd+Td)
        noise_level = sum(sum(db2pow(RDM([i-Tr-Gr:i-Gr-1 i+Gr:i+Gr+Tr-1],[j-Td-Gd:j-Gd-1 j+Gd:j+Gd+Td-1]))));
        thr = pow2db(noise_level/((2*(Td+Gd+1))^2 - (Gd*Gr) - 1));
        thr = thr + offset;
        
        Cut(i,j) = RDM(i,j);
        
        if Cut(i,j)<thr
            Cut(i,j)=0;
        else
            Cut(i,j)=1;
        end        
    end    
end


%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,Cut);
colorbar;


 
 