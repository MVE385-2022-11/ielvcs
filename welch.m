function w=welch(L)
    w = 1 - (2*(0:L-1)/(L-1) - 1).^2;
    w = w(:);
end
