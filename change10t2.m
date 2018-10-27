function code=change10t2(x)
N_row = length(x);
code = zeros(N_row,16); % 16-bit quant
code(x>=0,:) = de2bi(x(x>=0),16,'left-msb');
code(x<0,:) = de2bi(2^16-double(abs(x(x<0))),16,'left-msb');
end