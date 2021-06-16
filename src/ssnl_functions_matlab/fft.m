function field = fft(field)

% out = ifftshift(fft(fftshift(field))) ./ this.grid_nPts;
field = ifftshift(fft(fftshift(field)));

end