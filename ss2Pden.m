function [P, den] = ss2Pden(sys,tol)
% SS2PDELTA Converts a state-space system to a transfer function matrix
% written as Gp = 1/delta * P, where P is a polynomial matrix and delta is
% a polynomial.

% Get state-space function matrices
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

% Compute the transfer function matrix
s = zpk('s');
Gp = zpk(C * inv(s*eye(size(A)) - A) * B + D);
[r,m] = size(Gp);

% Place the zero at zero
for i = 1:r
    for j = 1:m
        for k = 1:numel(Gp.z{i,j})
            if abs(Gp.z{i,j}(k)) <= tol
                Gp.z{i,j}(k) = 0;
            end
        end
        for k = 1:numel(Gp.p{i,j})
            if abs(Gp.p{i,j}(k)) <= tol
                Gp.p{i,j}(k) = 0;
            end
        end
    end
end
Gp = minreal(Gp);

% Compute the denominator of the resolvent (LCM of all denominators)
lcmDen = [];
for i = 1:r
    for j = 1:m
        lcmDen = getLcm(lcmDen,Gp.p{i,j},tol);
    end
end
den = zpk(lcmDen,[],1);

% Compute the numerator matrix of the system
z    = cell(r,m);
z(:) = {[]};
p    = cell(r,m);
p(:) = {[]};
k    = zeros(r,m);
P = zpk(z, p, k);

for i = 1:r
    for j = 1:m
        this = Gp(i,j);
        z = getMissingZeros(this.p{1}, den.z{1}, tol);
        P.z{i,j} = [Gp.z{i,j} ; z];
%         P.k(i,j) = Gp.k(i,j) / prod(abs(z));
        P.k(i,j) = 1;
    end
end

% Evaluate the gain of the two transfer function at 1Hz to set the gain
old = evalfr(Gp,    1);
new = evalfr(P/den, 1);
for k = 1:numel(new)
    if new(k) ~= 0
        new(k) = old(k)/new(k);
    end
end
P.k = new;


end



function zC = getLcm(zA,zB,tol)
% GETLCM Return the zeros of the lowest common multipliers of two 
% polynomials whose zeros are zA and zB.

% Check if the zeros appears in both polynomial. If yes, remove it from A.
% If not, keep it, it will be added to the LCM.
for k = 1:numel(zB)
    idx = find(abs(zA - zB(k)) <= tol * abs(zA + zB(k)));
    if (~isempty(idx) && sum(idx) ~= 0)
        % Common zero
        zA(idx(1)) = [];
    end
end

zC = [zA;zB];

end



function zC = getGcd(zA,zB,tol)
% GETGCD Return the zeros of the greatest common divisors of two 
% polynomials whose zeros are zA and zB.

% Check if the zeros appears in both polynomial. If yes, this is a common
% zero.
zC = [];
for k = 1:numel(zB)
    idx = find(abs(zA - zB(k)) <= tol * abs(zA + zB(k)));
    if (~isempty(idx) && sum(idx) ~= 0)
        % Common zero
        zC = [zC; zB(k)];
    end
end

end



function z = getMissingZeros(this,den,tol)
% GETMISSINGPOLES Returns the zeros that appears in `den' but not in `this'.

% Compute the gcd of this and den.
gcd = getGcd(this,den,tol);

% Find the zeros that appear in den but not in the GCD and remove it.
z = den;
for k = 1:numel(gcd)
    idx = find(abs(z - gcd(k)) <= tol * abs(z + gcd(k)));
    z(idx(1)) = [];
end

end






