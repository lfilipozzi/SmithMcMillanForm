function [UL, UR, SMF] = smithMcMillanForm(P, den, sVar)
% SMITHMCMILLANFORM Compute the Smith-McMillan form and return the
% unimodular matrices UL and UR, such that SMF = UL*Gp*UR, where 
% Gp = P/den. All matrices are symbolics.
% 
% @author Louis Filipozzi

if nargin < 3
    sVar = symvar(P);
end

if isempty(sVar)
    sVar = symvar(tf2sym(tf('s')));
end
[UL, UR, SNF] = smithNormalForm(P, sVar);
SMF = simplify(SNF/den);

end



function [UL, UR, SNF] = smithNormalForm(P, sVar)
% SMITHNORMALFORM Compute the Smith Normal form of a symbolic matrix P and
% return the unimodular matrices so that SNF = UL*P*UR.

[r,m] = size(P);

if (m > r)
    % The matrix is wide.
    
    [ULt, URt, SNFt] = smithNormalForm(P.',sVar);
    SNF = SNFt.';
    UL  = URt.';
    UR  = ULt.';
else
    % The matrix is tall or square.
    
    % Compute the hermite form.
    [UH, H] = hermiteForm(P,sVar);
    
    % Extract square matrix.
    HS = H(1:m,1:m);
    
    % Compute the Smith normal form of the square matrix HS.
    if rank(HS) < m
        % The plant is singular.
        
        warning(strcat(...
            "The plant is singular. Use row and column manipulation ", ...
            "to diagonalize the plant. This may take some time."...
        ));
        
        % Use unimodular transformation to diagonalize the matrix.
        [ULD,URD,D] = unimodularDiagonalization(HS,sVar);
        d = rank(D);
        n = size(HS,1); % size of the square matrix HS
        % Compute the Smith Form of the nonsingular matrix.
        [ULNS, URNS, SNFNS] = smithForm(D(1:d,1:d));
        % Compute the Smith Form of the upper triangular square matrix HS.
        SNFS = blkdiag(SNFNS,zeros(n-d));
        ULS = simplify(blkdiag(ULNS,eye(n-d)) * ULD);
        URS = simplify(URD * blkdiag(URNS,eye(n-d)));
    else
        % The matrix is not singular.
        
        % Compute the Smith Form of the upper triangular square matrix HS.
        [ULS, URS, SNFS] = smithForm(HS,sVar);
    end
    
    % Return the Smith-McMillan form.
    SNF = sym(zeros(r,m));
    SNF(1:m,1:m) = SNFS;

    % Compute the unimodular matrices.
    UL = simplify(blkdiag(ULS,eye(r-m)) * UH);
    UR = URS;
end

end



function [UL,UR,D] = unimodularDiagonalization(P,sVar)
% UNIMODULARDIAGONALIZATION Compute left and right unimodular matrices to
% diagonalize a symbolic polynomial matrix. 
% 
% Return matrices UL, UR, and D such that UL*P*UR = D.
% 
% This function does not compute the Smith normal form as it only enforces
% monic polynomials on the diagonal with increasing degree. It does not
% enforce that each diagonal polynomial divides the following nonzero 
% polynomial.

[r,m] = size(P);

if (m < r)
    % The matrix is tall.
    
    [ULt,URt,Dt] = unimodularDiagonalization(P.',sVar);
    D = Dt.';
    UL  = URt.';
    UR  = ULt.';
else
    % The matrix is wide or square.
    
    % Create the augmented matrix.
    %        [ P | I ]   [  D  | UL ]
    % Paug = [---+---] ~ [-----+----]
    %        [ I |   ]   [  UR |    ]
    Paug = [...
        P      eye(r);
        eye(m) zeros(m,r);
    ];
    n = min(r,m);

    % Current row and column to cancel.
    jt = 1;

    % skipCol keep track of zero columns to ignore.
    skipCol = zeros(1,m);

    for t = 1:n
        % If the column is already only zero element AND if the column can 
        % be skipped (i.e., there are more remaining columns than entries 
        % to diagonalize: jt-t<col-n), then skip the column
        while isequal(Paug(1:r,jt),0*Paug(1:r,jt)) && jt < (m-n+t)
            skipCol(jt) = 1;
            jt = jt+1;
        end

        % As long as there is a non-zero 'off-diagonal' entry, use the
        % Euclidean division, row operation and permutaion to eliminate 
        % this entry.
        while not(isequal(Paug(t+1:r,jt), 0*Paug(t+1:r,jt)) && ...
                isequal(Paug(t,jt+1:m),0*Paug(t,jt+1:m)))
            if isa(Paug, 'nan')
                error("The matrix contains NaN");
            end

            % Step 1: Choose the pivot as the non-zero polynomial with the
            % lowest degree amongst the entries of the current row and 
            % column.
            
            lowest_deg = inf;
            for l = t:r
                deg_tmp = getDegree(Paug(l,jt),sVar);
                if deg_tmp < lowest_deg && ...
                        not(isequal(Paug(l,jt),sym(0)))
                    lowest_deg_row = l;
                    lowest_deg_col = jt;
                    lowest_deg = deg_tmp;
                end
            end
            for l = jt+1:m
                deg_tmp = getDegree(Paug(t,l),sVar);
                if deg_tmp < lowest_deg && ...
                        not(isequal(Paug(t,l),sym(0)))
                    lowest_deg_row = t;
                    lowest_deg_col = l;
                    lowest_deg = deg_tmp;
                end
            end

            % Step 2: Use row and column permutation to place the pivot.
            
            if lowest_deg_row ~= t
                % Row permutation
                Paug([t lowest_deg_row],:) = Paug([lowest_deg_row t],:);
            elseif lowest_deg_col ~= jt
                % Column permutation
                Paug(:,[jt lowest_deg_col]) = Paug(:,[lowest_deg_col jt]);
            end
            pivot = Paug(t,jt);

            % Step 3: Use Euclidian division and linear combination of rows
            % and columns to reduce the degree of all non-zero remaining 
            % polynomials of the current row and column.

            % Eliminate column entries.
            for l = t+1:r
                % Perform euclidean division.
                if ~isempty(symvar(pivot))
                    [Q,~] = quorem(Paug(l,jt),pivot,sVar);
                else
                    % Euclidian division when the pivot is a constant term.
                    Q = Paug(l,jt) / pivot;
                end
                % Use row operation: row <- row - quotient*pivot.
                Paug(l,:) = simplify(Paug(l,:) - Q*Paug(t,:));
            end
            % Eliminate row entries.
            for l = jt+1:m
                % Perform euclidean division.
                if ~isempty(symvar(pivot))
                    [Q,~] = quorem(Paug(t,l),pivot,sVar);
                else
                    % Euclidian division when the pivot is a constant term.
                    Q = Paug(t,l) / pivot;
                end
                % Use row operation: column <- column - quotient*pivot.
                Paug(:,l) = simplify(Paug(:,l) - Q*Paug(:,jt));
            end
        end

        % Make the pivot monic.
        pivot = Paug(t,jt);
        if ~isequal(pivot, sym(0))
            poly_tmp = sym2poly(pivot);
            coef = poly_tmp(1);
            Paug(:,jt) = Paug(:,jt)/coef;
            Paug = simplify(Paug);
        end

        % Move to the next row and column.
        jt = jt+1;
    end

    % Step 4: Rearrange the matrix to sort the diagonal entries by degree.
    % Rearrange to obtain a real diagonal in case a column has been 
    % skipped: if a column has been skiped, put it as the last column.
    
    I = find(skipCol);
    for l = length(I):-1:1
        Paug(:,1:m) = [Paug(:,1:I(l)-1) Paug(:,I(l)+1:m) Paug(:,I(l))];
    end
    % Sort the diagonal entries by decreasing degree
    deg_diag = zeros(1,n);
    for l = 1:n
        deg_diag(l) = getDegree(Paug(l,l),sVar);
    end
    % Set the zero polynomial at the end
    deg_diag(deg_diag == -inf) = inf;
    [~,I] = sort(deg_diag);
    Paug(:,1:m) = [Paug(:,I)  Paug(:,length(I)+1:m)];
    Paug(1:r,:) = [Paug(I,:); Paug(length(I)+1:r,:)];
    
    % Return the diagonal and unimodular matrices.
    D = Paug(1:r,1:m);
    UL = Paug(1:r, m+1:end);
    UR = Paug(r+1:end, 1:m);
end

end



function deg = getDegree(pol,sVar)
% GETDEGREE Return the degree of a polynomial. Return -inf if the
% polynomial is null.

pol = simplify(pol);

if isequal(pol, sym(0))
    deg = -inf;
else
    deg = polynomialDegree(pol, sVar);
end

% if isequal(pol,0*pol)
%     deg = -inf;
% else
%     deg = feval(symengine, 'degree', pol);
% end

end














