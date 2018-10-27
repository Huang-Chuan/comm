function sig=change2to10(code)
    N_row = size(code,1);
    sig = zeros(N_row,1);
    sig(code(:,1)==1) = bi2de(code(code(:,1)==1,:),'left-msb') - 2^16;
    sig(code(:,1)==0) = bi2de(code(code(:,1)==0,:),'left-msb');
end

