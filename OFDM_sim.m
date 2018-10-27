clear all;
close all;
B_subc = 15e3; %subcarrier
num_c = 4;  % 4 subcarrier

fc_0 = 15e3; % first carrier at 15KHz
fc_1 = 2*fc_0;
fc_2 = 3*fc_0;
fc_3 = 4*fc_0;


Ts = 1/fc_0;  %  OFDM symbol duration
fs = 1500e3;  %analog sampling rate 


t = 0:1/fs:(2*Ts-1/fs); % countinuous time, two OFDM symbol
N = length(t)/2;      % number of samples per OFDM symbol

QAM_symbol = [1+1i 1-1i; -1-1i -1-1i; -1+1i 1-1i; 1-1i -1+1i];

sig_c0 = [QAM_symbol(1,1)*ones(1,N),QAM_symbol(1,2)*ones(1,N)].*exp(j*2*pi*fc_0*t);
sig_c1 = [QAM_symbol(2,1)*ones(1,N),QAM_symbol(2,2)*ones(1,N)].*exp(j*2*pi*fc_1*t);
sig_c2 = [QAM_symbol(3,1)*ones(1,N),QAM_symbol(3,2)*ones(1,N)].*exp(j*2*pi*fc_2*t);
sig_c3 = [QAM_symbol(4,1)*ones(1,N),QAM_symbol(4,2)*ones(1,N)].*exp(j*2*pi*fc_3*t);

s = sig_c0+sig_c1+sig_c2+sig_c3;
plot(real(s));

tx = real(s);

%tx_I_c0 = tx.*cos(2*pi*fc_0*t);
%figure
%plot(tx_I_c0)
%u=[sum(tx_I_c0(1:100)),sum(tx_I_c0(101:200))]
 
%tx_Q_c0 = tx.*sin(2*pi*fc_0*t);
%figure
%plot(tx_Q_c0 )

%u=[sum(tx_Q_c0(1:100)),sum(tx_Q_c0(101:200))]
hold on;

OFDM_symbol = 4*reshape(ifft([QAM_symbol(4,:); QAM_symbol(1:3,:)]),[8,1]);
plot(1:25:7*25+1,real(OFDM_symbol),'o');
% t = 0:7;
% ts = linspace(0,7,200);
% [Ts,T] = ndgrid(ts,t);
% y = sinc(Ts - T)*(OFDM_symbol);
% figure
% plot(real(y))
%OFDM_symbol = reshape(ifft([QAM_symbol(4,:); QAM_symbol(1:3,:)]),[8,1]);
%u =s([1:25:1+7*25]);

% Fs = 44.1e3;
% filtertype = 'FIR';
% Fpass = 8e3;
% Fstop = 12e3;
% Rp = 0.1;
% Astop = 80;
% FIRLPF = dsp.LowpassFilter('SampleRate',fs, ...
%                              'FilterType',filtertype, ...
%                              'PassbandFrequency',30e3, ...
%                              'StopbandFrequency',Fstop, ...
%                              'PassbandRipple',Rp, ...
%                              'StopbandAttenuation',Astop);