clear all;
clc;
N_ofdm = 2048;
f_delta = 15e3;

N_subcarrier = 1320;
N_CP = 144;
Symbol_number = 14;
conv_poly = [23,35];
K = 5;
trel = poly2trellis(K,conv_poly);
tail = zeros(K-1,1);
mod_degree = 4;
code_rate = 0.5;
tb_len = 50;
bit_length = mod_degree*code_rate*N_subcarrier*Symbol_number-K+1;

SNR_db = 0:10:20;
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
%-----------------------audio recording-------------------------------


recObj = audiorecorder(8000, 16, 1);
disp('Press any key to start %g seconds of recording...\n');
pause;
disp('Start speaking.');
recordblocking(recObj, 2);
disp('End of Recording.');
audio_sound = getaudiodata(recObj,'int16');
%-----------------------int16 to binary data--------------------------

binary_data = change10t2(audio_sound);

sig = reshape(binary_data, size(binary_data,1)*size(binary_data,2),1);

%-----------------------Zero padding and packing----------------------
N_block = ceil(length(sig)/bit_length);
binary_pb = zeros(bit_length,N_block);
num_zeros = N_block*bit_length - length(sig);
transmit_sig = [sig; zeros(num_zeros,1)];
%---------------------------------------------------------------------

for loop_snr = 1 : length(SNR)
    err = 0;
    sigma = sqrt(1/SNR(loop_snr)/2);
    fprintf('Simulation initiates, SNR = %d \n ',SNR_db(loop_snr));

    for loop_block = 1 : N_block
        block_sig = transmit_sig((loop_block-1)*bit_length+1 : loop_block*bit_length);
        code = convenc([block_sig; tail],trel);
        interleaved = matintrlv(code,280,length(code)/280);
        symbol = qammod(interleaved,2^mod_degree,'UnitAveragePower',true,'InputType','bit');
        transmit_data = zeros(Symbol_number * (N_CP+N_ofdm), 1);
      
        for loop_symbol = 1 : Symbol_number
        	freq_domain = zeros(N_ofdm,1);
        	freq_domain((N_ofdm-N_subcarrier)/2+1:(N_ofdm-N_subcarrier)/2+N_subcarrier)= symbol((loop_symbol-1)*N_subcarrier+1:loop_symbol*N_subcarrier);
        	time_domain = ifft(freq_domain)*sqrt(N_ofdm);

            % add CP
            transmit_data((loop_symbol-1)*(N_CP+N_ofdm)+1:loop_symbol*(N_CP+N_ofdm)) = [time_domain(N_ofdm-N_CP+1:N_ofdm);time_domain];

            
        end        


        received_data = filter(channel, transmit_data);
        channel_coeff = channel.PathGains;

        received_data = zeros(length(transmit_data),1);
        for ii = 1 : length(received_data)
            for ll = 1 : length(tau)
                if(position(ll)< ii)
                received_data(ii) = received_data(ii) + transmit_data(ii - position(ll)) * channel_coeff(floor((ii-1)/(N_CP + N_ofdm)) * (N_CP + N_ofdm) + N_CP, ll);
                end
            end
        end


        noise = sigma*(randn(length(transmit_data),1) + j*randn(length(transmit_data),1));
 		received_data = received_data + noise;





        channel_delay_domain = zeros(N_ofdm,1);
 		for loop_symbol = 1 : Symbol_number

            channel_delay_domain(position+1) = channel_coeff((loop_symbol - 1)*(N_CP+N_ofdm) + N_CP,:);
        



 			de_cp = received_data((loop_symbol-1)*(N_CP+N_ofdm)+N_CP+1:loop_symbol*(N_CP+N_ofdm));
            channel_freq_domain = fft(channel_delay_domain);
 			fft_data = (fft(de_cp).*conj(channel_freq_domain))./(channel_freq_domain.*conj(channel_freq_domain))/sqrt(N_ofdm);        % this is frequency domain




 			demapp_data((loop_symbol-1)*N_subcarrier+1:loop_symbol*N_subcarrier) = fft_data((N_ofdm-N_subcarrier)/2+1:(N_ofdm-N_subcarrier)/2+N_subcarrier);

            noise_var((loop_symbol - 1) * N_subcarrier + 1 : loop_symbol * N_subcarrier) = sigma^2 ./ abs(channel_freq_domain((N_ofdm - N_subcarrier)/2+1:...
                   (N_ofdm - N_subcarrier)/2+ N_subcarrier )).^2;
        end

        fprintf('decoding......\n');
        for nsym = 1 : length(demapp_data)

    	data_demodulated(:,nsym) = qamdemod(demapp_data(nsym),2^mod_degree,'UnitAveragePower',true,'OutputType','approxllr','NoiseVariance',noise_var(nsym));
        end

    	temp1 = size(data_demodulated);

    	bit_stream = reshape(data_demodulated,temp1(1)*temp1(2),1);
        deinterleaved = matdeintrlv(bit_stream, 280, length(code)/280);


    	decision = vitdec(deinterleaved, trel, tb_len, 'term', 'unquant');
    	decision = decision(1:bit_length);

        binary_pb(:,loop_block) = decision; %store 
        
%         err = err + sum(decision ~= block_sig);
%         fprintf('block err = %d \n ',err);
%         
%         err_ = err + sum(binary_pb(loop_block,:)~= block_sig');
%         fprintf('block err_ = %d \n ',err);
    	% if(sum(decision~=sig)~=0)
    	% 	err_blk = err_blk + 1;
    	% end

    	% if(err_blk>=10)
    	% 	break;
    	% end
    end
    % ber(loop_snr) = err/(bit_length*loop_block);
    % bler(loop_snr) = err_blk/loop_block;
    % fprintf('Total bit error rate = %f \n',ber(loop_snr));
    
    
    decoded = reshape(binary_pb(1:length(sig))',[], 16);
    receive_audio = int16(change2to10(decoded));
    player = audioplayer(receive_audio,8000);
    play(player);
end


   

% figure;
% plot(real(fft_data),imag(fft_data),'r.')
% hold on
% plot(real(symbol),imag(symbol),'ko')
% axis([-1 1 -1 1])


% figure;
% semilogy(SNR_db,ber,'-^')
% grid on
% xlabel('SNR(dB)')
% ylabel('BER')        