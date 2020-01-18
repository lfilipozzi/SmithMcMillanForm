clear; clc;
close all;

%% Select the plant to test the function
% Choose plant:
% 1: Random plant
% 2: Plant from the example of chapter 13
% 3: Anti-roll bar system used by Sergio
% 4: Half-car model with load transfer and normal forces as inputs
select_plant = 3;

switch select_plant
    case 1
        % Random plant
        P = [ ...
            tf([1 1 3 2 4],1) tf([1 2 3 5],1) tf([1 3],1);
            tf([2 1 4 2 3],1) tf([2 2 2 3],1) tf([2 3],1);
        ];

        den = tf([1 5 4 4 3 5 3],1);

        Gp = P/den;
        
    case 2
        % Example from chapter 13
        P = 100 * [ ...
            tf([1 0],1) tf([1 0],1) tf([0 1],1);
            tf([0 1],1) tf([0 1],1) tf([0 1],1);
        ];

        den = tf([1 1],1);

        Gp= P / den;
        
    case 3
        % Anti-roll bar
        ms = 1632.93/2;
        mu = 163.29/2;
        f = 1;
        wn = 2*pi*f;
        zeta = .5;
        k = .5*ms*wn*wn;
        kt = 10*k;
        b = .5*zeta*wn*ms;
        bt = b/100;
        w = 1.825;
        J = 226.5241; 
        h = .762; %m
        l = .45; %m
        ktau = 85*180/pi; %Nm/rad
        btd = b/100;
        % Plant description
        A = [0 0 0 0 0 1/J 0 0 0 0;
            0 0 -1/mu 0 0 0 0 0 0 0;
            0 kt -(b+bt)/mu -k b/ms -b*w*.5/J 0 0 0 ktau/l;
            0 0 1/mu 0 -1/ms w*.5/J 0 0 0 0;
            0 0 b/mu k -2*b/ms 0 k b/mu 0 0;
            0 0 -b*w*.5/mu -k*w*.5 0 -b*w*w*.5/J k*w*.5 b*w*.5/mu 0 0;
            0 0 0 0 -1/ms -w*.5/J 0 1/mu 0 0;
            0 0 0 0 b/ms b*w*.5/J -k -(b+bt)/mu kt -ktau/l;
            0 0 0 0 0 0 0 -1/mu 0 0;
            0 0 -1/(l*mu) 0 0 0 0 1/(l*mu) 0 0];
        B = [0 0 0;0 0 0;-1 0 0;0 0 0;1 1 0;-w/2 w/2 0;0 0 0;0 -1 0;0 0 0;0 0 1];   %Fdl;Fdr;wm
        Bd= [0 0;1 0;btd 0;0 0;0 0;0 ms*h;0 0;btd 0;1 0;0 0]; %dxi;ay
        alpha = .99999;
        Kaa = [.5 .5 0;alpha/2 -alpha/2 (1 - alpha)];
        Ba = pinv(Kaa);
        %kh = .5;
        %kr2 = .001;
        %kr1 = kr2*5500;
        %Ba= [kh kr1;kh -kr1;0 kr2];
        C = [0 0 0 0 1/ms 0 0 0 0 0;    %dxs
            1 0 0 0 0 0 0 0 0 0];   %phi
        D = [0 0 0;0 0 0];
        Da = [0 0;0 0];
        Dd = [0 0;0 0];
    %     sys = ss(A,B,C,D);
    %     sys_d = ss(A,Bd,C,Dd);
    %     % Tranfser Function Description
    %     Gd = minreal(tf(sys_d),1e-03);
    %     G = minreal(tf(sys),1e-03);

        s = tf('s');
        resolvent = inv(s * eye(size(A)) - A);
    %     for k = 2:10
    %         resolvent.num{1,k} = resolvent.num{1,k}(1:end-1);
    %         resolvent.den{1,k} = resolvent.den{1,k}(1:end-1);
    %     end
        den = tf(resolvent.den{1,2},1);
    %     for k = 1:10
    %         resolvent.num{k,1} = den.num{1}(1:end-1);
    %         resolvent.den{k,1} = den.num{1};
    %     end
        P = C * tf(resolvent.num,1) * B + D * den;


        % Place the zero at zero
        tol = 1e-5;
        for i = 1:size(P,1)
            for j = 1:size(P,2)
                for k = flip(1:numel(P.num{i,j}))
                    if abs(P.num{i,j}(k)) <= tol
                        P.num{i,j}(k) = 0;
                    else
                        break;
                    end
                end
            end
        end

        for k = flip(1:numel(den.num{1}))
            if abs(den.num{1}(k)) <= tol
                den.num{1}(k) = 0;
            else
                break;
            end
        end

        % Find common term in all transfer functions
        den = zpk(den);
        P = zpk(P);
        i = 1;
        while i < numel(den.z{1})
            % Check if the zero of the denominator is a zero of all TF
            isZInTFM = true;
            for j = 1:numel(P.z)
                if ~ismember(den.z{1}(i), P.z{j})
                    isZInTFM = false;
                    break;
                end
            end

            % The zero of the denominator can be reomved from the TFM.
            if isZInTFM
                % Message
                disp(strcat(num2str(den.z{1}(i)), " is a canceled zero-pole."));
                % Remove from the numerator matrix
                for j = 1:numel(P.z)
                    idx = find(den.z{1}(i) == P.z{j});
                    % Remove only one zero (in case of multiple zero)
                    P.z{j}(idx(1)) = [];
                end
                % Remove from the denominator
                den.z{1}(i) = [];
                % Reindex counter
                i = i-1;
            end
            % Increment counter
            i = i+1;
        end

        Gp = P/den;
    case 4
        % Half-car model
        ms = 980;
        mu = 50;
        ks = 25000;
        bs = 1000;

        A = [...
            -bs/ms,  bs/mu,  ks;
             bs/ms, -bs/mu, -ks;
             -1/ms,   1/mu,   0;
        ];

        B = [...
            0, -1;
            1,  0;
            0,  0;
        ];

        C = [...
            -bs/ms^2, bs/(ms*mu), ks/ms;
                   0,          0,     1;
               -1/ms,       1/mu,     0;
        ];

        D = [...
            0, -1/ms;
            0,     0;
            0,     0; 
        ];

        s = tf('s');
        resolvent = inv(s * eye(size(A)) - A);
        den = tf(resolvent.den{1,1},1);
        P = C * tf(resolvent.num,1) * B + D * den;


        % Place the zero at zero
        tol = 1e-8;
        for i = 1:size(P,1)
            for j = 1:size(P,2)
                for k = flip(1:numel(P.num{i,j}))
                    if abs(P.num{i,j}(k)) <= tol
                        P.num{i,j}(k) = 0;
                    else
                        break;
                    end
                end
            end
        end

        for k = flip(1:numel(den.num{1}))
            if abs(den.num{1}(k)) <= tol
                den.num{1}(k) = 0;
            else
                break;
            end
        end

        % Find common term in all transfer functions
        den = zpk(den);
        P = zpk(P);
        i = 1;
        while i < numel(den.z{1})
            % Check if the zero of the denominator is a zero of all TF
            isZInTFM = true;
            for j = 1:numel(P.z)
                if ~ismember(den.z{1}(i), P.z{j})
                    isZInTFM = false;
                    break;
                end
            end

            % The zero of the denominator can be reomved from the TFM.
            if isZInTFM
                % Message
                disp(strcat(num2str(den.z{1}(i)), " is a canceled zero-pole."));
                % Remove from the numerator matrix
                for j = 1:numel(P.z)
                    idx = find(den.z{1}(i) == P.z{j});
                    % Remove only one zero (in case of multiple zero)
                    P.z{j}(idx(1)) = [];
                end
                % Remove from the denominator
                den.z{1}(i) = [];
                % Reindex counter
                i = i-1;
            end
            % Increment counter
            i = i+1;
        end

        Gp = P/den;
        
    case 5
        % Rank deficient plant
        P = 100 * [ ...
            tf([1 2],1) tf([1 2],1) tf([1 2],1);
            tf([1 2],1) tf([1 2],1) tf([1 2],1);
        ];

        den = tf([1 1],1);

        Gp= P / den;

    case 6
        % Rank deficient plant
        P = 100 * [ ...
            tf([1 2],1) tf([1 2],1) tf([1 2],1) tf([1 4],1) tf([1 2],1) tf([1 2],1);
            tf([1 2],1) tf([1 2],1) tf([1 2],1) tf([1 2],1) tf([1 2],1) tf([1 2],1);
            tf([1 2],1) tf([1 2],1) tf([1 2],1) tf([1 5],1) tf([1 2],1) tf([1 2],1);
            tf([1 2],1) tf([1 2],1) tf([1 2],1) tf([1 2],1) tf([1 2],1) tf([1 2],1);
        ];

        den = tf([1 1],1);

        Gp= P / den;

    case 7
        % Rank deficient plant
        P = 100 * [ ...
            tf([0 2],1) tf([1 2],1) tf([1 2],1) tf([1 4],1) tf([1 2],1) tf([1 2],1);
            tf([0 2],1) tf([1 2],1) tf([1 2],1) tf([1 2],1) tf([1 2],1) tf([1 2],1);
            tf([0 2],1) tf([1 2],1) tf([1 2],1) tf([1 5],1) tf([1 2],1) tf([1 2],1);
            tf([0 2],1) tf([1 2],1) tf([1 2],1) tf([1 2],1) tf([1 2],1) tf([1 2],1);
        ];

        den = tf([1 1],1);

        Gp= P / den;

    case 8
        % Rank deficient plant
        s = tf('s');
        P = 100 * [ ...
            s^3 s^2 s^1;
            0   0   s^2;
            0   0   s^3;
        ];

        den = tf([1 1],1);

        Gp= P / den;

    case 9
        % Rank deficient plant
        s = tf('s');
        P = 100 * [ ...
            0 0 s;
            0 0 0;
            0 0 0;
        ];

        den = tf([1 1],1);

        Gp= P / den;

    case 10
        % Rank deficient plant
        s = tf('s');
        P = 100 * [ ...
            0   s^3 s^2 s^1;
            0   0   0   s^2;
            0   0   0   s^3;
            0   0   0   s^4
        ];

        den = tf([1 1],1);

        Gp= P / den;
        
    otherwise
        error('No plant');
end

%% Test function
Psym   = tf2sym(P);
densym = tf2sym(den);

% Psym   = vpa(Psym,3);
% densym = vpa(densym,3);

[ULsym, URsym, MPsym] = smithMcMillanForm(Psym, densym);

simplify(ULsym * Psym * URsym / densym - MPsym)
vpa(det(ULsym),2)
vpa(det(URsym),2)

MP = zpk(sym2tf(MPsym));
UL = zpk(sym2tf(ULsym));
UR = zpk(sym2tf(URsym));


