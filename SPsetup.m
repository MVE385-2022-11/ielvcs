classdef SPsetup
    properties
        L; % MTI filter tap count
        M; % Doppler freq resolution

        h; % FIR
        H; % FIT convolution matrix
        W; % Windowing weights
        B; % Freq base matrix BN
        C; % Unequalized SP
        F; % Equalized SP
    end
    methods
        % Constructor
        function obj = SPsetup(smod, varargin)
            % varargs
            if nargin == 1
                obj.h = [-1 2 -1]/sqrt(6);
            else
                obj.h = varargin{1};
            end

            obj.L = length(obj.h);
            obj.M = smod.N - obj.L + 1; % == size(H,1)

            if nargin <= 2
                win = hamming(obj.M);
            else
                win = varargin{2};
            end

            %% MTI
            % FIR choice
            %h = [1];
            %h = [-1 1]/sqrt(2);
            %h = [-1 2 -1]/sqrt(6);
            %h = normalize([-1 16 -30 16 -1], 'norm', 2);  % Experiment
            %h = normalize([1 -4 6 -4 1], 'norm', 2);  % Experiment
            obj.H = convmtx(fliplr(obj.h),obj.M);
            % equivalent to
            % conv(y,h, 'valid');
            
            %% WDFT
            %ham = 25/46 - (1 - 25/46) * cos(2*pi*(0:M-1)/M);
            %W = diag(ham);
            %W = diag(hamming(M)); % slightly different
            %W = eye(M);
            %W = diag(hann(M));
            %W = diag(kaiser(M,20));
            obj.W = diag(win);
            
            %B = b(M).^repmat(fp * (0:M - 1)/M, M, 1);
            obj.B = dftmtx(obj.M)';
            %fy = B' * W * y;
            % equivalent to
            %fy = fft(ham' .* y);
            
            %% EQ
            obj.C = obj.B' * obj.W * obj.H;
            e = vecnorm(obj.C,2,2);
            obj.F = obj.C ./ e;
        end
    end
end
