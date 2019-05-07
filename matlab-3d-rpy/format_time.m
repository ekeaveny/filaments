function [elapsed_time_hms] = format_time(elapsed_time)
% FORMAT_TIME  Formats a time in seconds as something more readable and
%              (mostly) fixed width.
%
%   FORMAT_TIME(20)     returns '   20.0s'
%   FORMAT_TIME(200)    returns '    3:20'
%   FORMAT_TIME(20000)  returns ' 5:33:20'
%   FORMAT_TIME(200000) returns '2d 5:33:20'

    if elapsed_time >= 86400
        tr2b_m = floor(elapsed_time/60);
        tr2b_s = mod(elapsed_time,60);
        tr2b_h = floor(tr2b_m/60);
        tr2b_m = mod(tr2b_m,60);
        tr2b_d = floor(tr2b_h/24);
        tr2b_h = mod(tr2b_h,25);
        elapsed_time_hms = sprintf('%dd%2d:%02d:%02.0f', tr2b_d, tr2b_h, tr2b_m, tr2b_s);
    elseif elapsed_time >= 3600
        tr2b_m = floor(elapsed_time/60);
        tr2b_s = mod(elapsed_time,60);
        tr2b_h = floor(tr2b_m/60);
        tr2b_m = mod(tr2b_m,60);
        elapsed_time_hms = sprintf('%2d:%02d:%02.0f', tr2b_h, tr2b_m, tr2b_s);
    elseif elapsed_time >= 60
        tr2b_m = floor(elapsed_time/60);
        tr2b_s = mod(elapsed_time,60);
        tr2b_h = floor(tr2b_m/60);
        tr2b_m = mod(tr2b_m,60);
        elapsed_time_hms = sprintf('   %2d:%02.0f',tr2b_m, tr2b_s);
    else
        elapsed_time_hms = sprintf('%7.1fs',elapsed_time);
    end
end