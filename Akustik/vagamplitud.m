
r_start = 0.1;
r = r_start:0.01:5;
amp = 4*pi;
amps = zeros(length(r),1);

for i = 1:length(r)
    amps(i) = amplitude(r(i), amp);
end

plot(r, amps, 'k')
xlabel('Distance r')
ylabel('Sound Pressure')
xline(1, '--k');
title('Sound Pressure As a Function of Distance')
legend('Sound Pressure','Standard Distance of Reference')


function a = amplitude(r, amp)
    a = amp/(4*pi*r);
end