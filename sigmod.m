classdef sigmod
    properties(Constant)
        % speed of light
        c = physconst('lightspeed');
    end
    properties
        % Radar specs
        N; % # of pulses
        f0; %Hz, carrier frequency
        fp; %Hz, pulse repetition frequency

        % Useful constants
        vp;
        b0;
    end
    methods
        % Constructor
        function obj = sigmod(num_pulses, carrier_freq, prf)
            % Radar specs
            obj.N = num_pulses; % # of pulses
            obj.f0 = carrier_freq; %Hz, carrier frequency
            obj.fp = prf; %Hz, pulse repetition frequency
            obj.vp = obj.f2v(obj.fp);
            obj.b0 = obj.b(obj.N);
        end
        % Doppler shift
        function fs = v2f(obj, vs)
            fs = mod(2 * vs/obj.c * obj.f0, obj.fp);
        end
        function vs = f2v(obj, fs)
            vs = fs * obj.c / (2 * obj.f0);
        end

        % steering vector base (take elementwise power of freq)
        function b0 = b(obj, M)
            b0 = exp(1i*2*pi*1/obj.fp * (0:M-1)).';
        end

        % create a signal from amplitudes and velocities
        function y = s(obj,as,vs)
            y = (obj.b0).^repmat(obj.v2f(vs), obj.N, 1) * as(:);
        end

        % sample complex circularly symmetric Gaussian
        function n = randcn(obj, sd, varargin)
            switch nargin
                case 2
                    sz = [obj.N 1];
                case 3
                    sz = [varargin(1) 1];
                otherwise
                    sz = cell2mat(varargin);
            end
            n = sd/sqrt(2) * (randn(sz) + 1i * randn(sz));
        end

        % clutter simulator
        function c=clutter(obj,NC,varargin)
            % Idea 1:
            % Generate random amplitudes and velocities
            %as = 1e5 .* (1 + 0.5 .* rand(1,NC))
            %vs = 2.*randn(1,NC)
            
            optargin = {1, 1e3, 4};
            optargin(1:(nargin - 2)) = varargin;
            K = optargin{1}; 
            max_a = optargin{2};
            max_v = optargin{3};

            %as = exprnd(1,1,NC);
            as = unifrnd(-max_a,max_a, K,NC);
            vs = unifrnd(-max_v,max_v, K,NC);
            
            c = zeros(obj.N,K);
            for i = 1:K
                c(:,i) = obj.s(as(i,:), vs(i,:));
                c(:,i) = c(:,i) ./ max(abs(c(:,i))) * max_a;
            end
            
            
            % Idea 2:
            %lowpass(sum(1e6/sqrt(2) * (randn(N,NC) + 1i * randn(N,NC)), 2),fp/(1e10*N),fp)
            %....
        end

    end
end
