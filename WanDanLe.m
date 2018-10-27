clear all;
clc;

%----------------------- OFDM modulation ------------------------
N_ofdm = 2048;                % IFFT 
f_delta = 15e3;               % subcarrier spacing(NOT used in simulation)
N_block = 1e3;                % maximum transmitting block
N_subcarrier = 1320;          % Number of subcarrier used in transmitting data
N_CP = 144;                   % number of samples as cyclic prefix 
Symbol_number = 14;           % number of 

%------------------------ convolution code------------------------
conv_poly = [23,35];  % g1 =[ 1 0 0 1 1] g2 = [1 1 1 0 1]
K = 5;                % constraint length
trel = poly2trellis(K,conv_poly);
tail = zeros(K-1,1);


%------------------------  frame ---------------------------------

mod_degree = 4;               % 16 QAM
code_rate = 0.5;              % 1/2 code rate
tb_len = 50;                  % trace back length used in viterbi decoding
bit_length = mod_degree*code_rate*N_subcarrier*Symbol_number-K+1;

SNR_db = 10:2:18;
SNR = 10.^(SNR_db/10);

ber = zeros(1,length(SNR));
bler = zeros(1,length(SNR));


%------------------Channel setting------------------------------------
ts = 1/30.72e6;
fd = 3;
tau = [0 310 710 1090 1730 2510] * 1e-9;
position = floor(tau/ts);
pd = [0 -1 -9 -10 -15 -20];
channel = rayleighchan(ts,fd,tau,pd);
channel.ResetBeforeFiltering = 1;
channel.StorePathGains = 1;
%---------------------------------------------------------------------

for loop_snr = 1 : length(SNR)
    err = 0;
    err_blk = 0;
    sigma = sqrt(1/SNR(loop_snr)/2);
    fprintf('Simulation initiates, SNR = %d \n ',SNR_db(loop_snr));

    for loop_block = 1 : N_block
        sig = round(rand(bit_length,1));
        code = convenc([sig; tail],trel);                    % convolution coding 
        interleaved = matintrlv(code,280,length(code)/280);  % interleaving 
        symbol = qammod(interleaved,2^mod_degree,'UnitAveragePower',true,'InputType','bit');
        transmit_data = zeros(Symbol_number * (N_CP+N_ofdm), 1);  % reserve space for 14 OFDM symbols
      
        for loop_symbol = 1 : Symbol_number
        	freq_domain = zeros(N_ofdm,1);
        	freq_domain((N_ofdm-N_subcarrier)/2+1:(N_ofdm-N_subcarrier)/2+N_subcarrier)= symbol((loop_symbol-1)*N_subcarrier+1:loop_symbol*N_subcarrier);   
            % encode QAM symbols into 1320 subcarriers which sit i  
        	time_domain = ifft(freq_domain)*sqrt(N_ofdm);    % convert into time domain
            % ifft(X) = 1/N_ofdm * (.....)
            % scaling term sqrt(N_ofdm) is to make sure time domain signal energy equal to QAM signal energy


            % add CP
            transmit_data((loop_symbol-1)*(N_CP+N_ofdm)+1:loop_symbol*(N_CP+N_ofdm)) = [time_domain(N_ofdm-N_CP+1:N_ofdm);time_domain];

            
        end        


        received_data = filter(channel, transmit_data);        % this is just to generate channel coefficients, filtering results NOT used
        channel_coeff = channel.PathGains;


        received_data = zeros(length(transmit_data),1);        % simulate channel output
        for ii = 1 : length(received_data)
            for ll = 1 : length(tau)
                if(position(ll)< ii)
                received_data(ii) = received_data(ii) + transmit_data(ii - position(ll)) * channel_coeff(floor((ii-1)/(N_CP + N_ofdm)) * (N_CP + N_ofdm) + N_CP, ll);
                % we assume channel DOES NOT change during one OFDM symbol period
                % linear convolution
                end
            end
        end


        noise = sigma*(randn(length(transmit_data),1) + j*randn(length(transmit_data),1));    % Gaussian noise 
 		received_data = received_data + noise;                                                % received signal 





        channel_delay_domain = zeros(N_ofdm,1);
 		for loop_symbol = 1 : Symbol_number

            channel_delay_domain(position+1) = channel_coeff((loop_symbol - 1)*(N_CP+N_ofdm) + N_CP,:);
    

 			de_cp = received_data((loop_symbol-1)*(N_CP+N_ofdm)+N_CP+1:loop_symbol*(N_CP+N_ofdm));
            channel_freq_domain = fft(channel_delay_domain);
 			fft_data = (fft(de_cp).*conj(channel_freq_domain))./(channel_freq_domain.*conj(channel_freq_domain))/sqrt(N_ofdm);        % this is frequency domain
            % Equalization in frequency domain 
            % scaling factor  sqrt(N_ofdm) again is to make sure equal energy 



 			demapp_data((loop_symbol-1)*N_subcarrier+1:loop_symbol*N_subcarrier) = fft_data((N_ofdm-N_subcarrier)/2+1:(N_ofdm-N_subcarrier)/2+N_subcarrier);

            noise_var((loop_symbol - 1) * N_subcarrier + 1 : loop_symbol * N_subcarrier) = sigma^2 ./ abs(channel_freq_domain((N_ofdm - N_subcarrier)/2+1:...
                   (N_ofdm - N_subcarrier)/2+ N_subcarrier )).^2;
            %noise variance is calculated at each carrier frequency used for qam demodulation 
            %white noise through linear filter
            %energy spectrum is scaled by |H(f)|^2
        end

        fprintf('decoding......\n');
        for nsym = 1 : length(demapp_data)

    	data_demodulated(:,nsym) = qamdemod(demapp_data(nsym),2^mod_degree,'UnitAveragePower',true,'OutputType','approxllr','NoiseVariance',noise_var(nsym));
        end

    	temp1 = size(data_demodulated);

    	bit_stream = reshape(data_demodulated,temp1(1)*temp1(2),1);
        deinterleaved = matdeintrlv(bit_stream, 280, length(code)/280);       % deinterleaving
    	decision = vitdec(deinterleaved, trel, tb_len, 'term', 'unquant');    % viterbi decoding, soft decision

    	decision = decision(1:bit_length);
    	err = err + sum(decision ~= sig);

    	if(sum(decision~=sig)~=0)
    		err_blk = err_blk + 1;
    	end

    	if(err_blk>=10)
    		break;
    	end
    end
    ber(loop_snr) = err/(bit_length*loop_block);
    bler(loop_snr) = err_blk/loop_block;
    fprintf('Total bit error rate = %f \n',ber(loop_snr));
    
    
    figure;
    plot(real(demapp_data),imag(demapp_data),'r.')
    hold on
    plot(real(symbol),imag(symbol),'ko')
    axis([-2 2 -2 2])

end



figure;
semilogy(SNR_db,ber,'-^')
grid on
xlabel('SNR(dB)')
ylabel('BER')        