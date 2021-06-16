function field = ifft(field)

% out = ifftshift(ifft(fftshift(field))) .* this.grid_nPts;
field = ifftshift(ifft(fftshift(field)));

end