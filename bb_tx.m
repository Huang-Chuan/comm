function r =  bb_tx(data,Ns,pulse)
% [r,data] =  bb_tx(Nbits,SIR,fj,SNR,Ns,pulse)
% Baseband binary signal source including noise and jamming
%
% ----------------- Inputs ---------------------------------
% Nbits = number of bits to simulate
%   SIR = signal-to-interference (jammer) ratio in dB
%    fj = bit rate normalized jammer frequency
%   SNR = signal-to-noise ratio in dB
%    Ns = number of samples per bit in waveform simulation
% pulse = transmitted pulse type, RECT = rectangular, 
%         SRC = square-root raised cosine with alpha = 0.5
%
% ----------------- Outputs ----------------------------------
%    r = transmitted signal as a data structure containing variations
%        of signal (S), jamming (J), and noise (N), in particular
%        r.SJN = S + J + N
%        r.SJ = S + J (signal plus jamming only)
%        r.N = N (noise only)
%        r.SNR = SNR to use later inside the receiver bb_rx
%        r.NS = NS to use later inside the receiver bb_rx
%        r.pulse = pulse to use later inside the receiver bb_rx
% data = data bits transmitted as +1/-1 bit values; can be used for
%        error checking in Monte-Carlo type simulation
% -------------------------------------------------------------
%
% Mark Wickert, November 2006


%data = sign(rand(1,Nbits)-0.5);
x = [data; zeros(Ns-1,length(data))];
x = reshape(x,1,Ns*length(data));

switch lower(pulse) % not case sensitive
    case 'src'
        h_src = sqrt_rc_imp(Ns,0.5,6);
        x = filter(h_src,1,x);
    case 'rect'
        h_rect = ones(1,Ns);
        x = filter(h_rect,1,x);
    otherwise
        error('pulse must be SRC or RECT')    
end


% Create a single tone jammer
% n = 0:length(x)-1;
% x_jammer = cos(2*pi*fj/Ns*n + 2*pi*rand(1,1));
% x_jammer = sqrt(Px*10^(-SIR/10))*x_jammer;
% Create white Gaussian noise
% w = sqrt(Ns*Px*10^(-SNR/10))*randn(size(x));
% Form output signal combinations in data structure variable r
r.SJN = x;  %+ x_jammer + w; % Complete signal
%r.SJ = x + x_jammer; % Signal + Jammer reference signal
%r.N = w;              % Noise alone reference signal
%r.SNR = SNR;
r.Ns = Ns;
r.pulse = pulse;
