
r_start = 0.1;
r = r_start:0.01:7;
amp = 4*pi;
amps = zeros(length(r),1);

for i = 1:length(r)
    amps(i) = amplitude(r(i), amp);
end

plot(r, amps, 'k')
xlabel('Distance r')
ylabel('Energy, Sound Pressure')
title('Amplitude of sound pressure as a function of distance from a point source')


function a = amplitude(r, amp)
    a = amp/(4*pi*r);
end