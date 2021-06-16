function fwhm = FWHM(lst)

if ~isreal(lst)
    lst = abs(lst).^2;
end

halfMax = max(lst)/2;
logicMat = find(lst > halfMax);
fwhm = length(logicMat(1):logicMat(end));


end