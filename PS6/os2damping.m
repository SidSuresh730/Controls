function f = os2damping(os)
    f = sqrt(log(os/100)^2/(pi^2 + log(os/100)^2));
end