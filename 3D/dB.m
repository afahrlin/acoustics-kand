function v = dB(u)
    p0 = 20*10^(-6);
    v = 20*log10(abs(u)./p0);
end