function [output] = SOFdis(Hs)
    [num, den] = tfdata(Hs, 'v');
    if den(1) ~= 1
            num = num/den(1);
            den = den/den(1);
    end
    Hs = tf(num, den);
    csys = canon(Hs, 'companion');
    A = rot90(csys.A,2);
    B = [];
    for i = 1:length(num)-1
        B(i) = num(i+1)-den(i+1)*num(1);
    end
    C = csys.B.';
    D = csys.D;
    
    output = ss(A,B.',C,D);
end