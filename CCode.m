function [output1] = CCode( m,g1,g2 ) 
    msg = [m ];
    n = length(msg);
    k = length(g1)-1;
    state = zeros(1,k);
    for i = 1:1:n
        m1(i) = rem(g1*[msg(i) state]',2); 
        m2(i) = rem(g2*[msg(i) state]',2);
        state = [msg(i) state(1:k-1)];
    end
   
    for i = 1:1:n
        output1(2*i-1) = rem(m1(i),2);
        output1(2*i) = rem(m2(i),2);
    end     
    
    output1 = output1.'; 
end

