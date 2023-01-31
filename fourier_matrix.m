function B=fourier_matrix(N, Nprime)
    f1 = (0:N-1).';
    f2 = (0:Nprime - 1).'/Nprime;
    % the fp frequency cancelling between the two.

    B = exp(1i*2*pi .* (f1*f2.'));
end
