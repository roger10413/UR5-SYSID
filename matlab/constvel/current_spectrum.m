function [order_axis, mag, f_motor] = current_spectrum(p, fs, n_gear)
    x = p.i - mean(p.i);
    n = numel(x);
    w = 0.5 - 0.5*cos(2*pi*(0:n-1)'/(n-1));  % Hanning window，Octave/MATLAB通用寫法
    xw = x .* w;
    Y = fft(xw);
    half = floor(n/2) + 1;
    f = (0:half-1)' * fs / n;
    mag = abs(Y(1:half)) * 2 / sum(w);
    f_motor = p.level * n_gear / (2*pi);   % 馬達端每秒圈數
    order_axis = f / f_motor;
end
