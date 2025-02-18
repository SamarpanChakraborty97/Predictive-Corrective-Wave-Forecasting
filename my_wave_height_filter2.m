function [wave_height,wave_idx]=my_wave_height_filter2(hz)
len_hz = length(hz);
zero_crossing(1) = 1;
for i=2:len_hz-1
    if (sign(hz(i+1) - hz(i))~=(sign(hz(i)-hz(i-1))))
        zero_crossing(end+1) = i;
    end
end

wave_idx = zero_crossing;
for i=1:length(wave_idx)
    wave_height(i) = hz(wave_idx(i));
end
end
        