function [Hw,num,den] = myDiscretize(Hs, Ts, approx) % default TUSTIN
    %s = tf('s');
    %z = tf('z', Ts);
    csys = SOFdis(Hs);
    if strcmp(approx, 'backward')
        Aw = inv(eye(length(csys.A)) - Ts*csys.A);
        Bw = Aw*csys.B*Ts;
    elseif strcmp(approx, 'forward')
        Aw = eye(length(csys.A)) + Ts*csys.A;
        Bw = csys.B*Ts;
    else % TUSTIN %
        Aw = inv(eye(length(csys.A)) - Ts/2*csys.A)*(eye(length(csys.A)) + Ts/2*csys.A);
        Bw = inv(eye(length(csys.A)) - Ts/2*csys.A)*csys.B*Ts/2;
    end
    dsys = ss(Aw, Bw, csys.C, csys.D, Ts);
    [b, a] = ss2tf(dsys.A, dsys.B, dsys.C, dsys.D);
    Beta = [];
    Alpha = [a(1)];
    for i = 1:length(b)-1
        Beta = [Beta, b(i)+b(i+1)-dsys.D*a(i+1)]; % DIFF FROM NOTES
        Alpha = [Alpha, a(i+1)];
    end
    Beta = [Beta, b(length(b))]; % DIFF FROM NOTES
    Hw = tf(Beta, Alpha, Ts);
end

    