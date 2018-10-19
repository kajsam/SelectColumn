function [w, h, Z] = select_column_all(X,Z,mask)

% Kajsa Mollersen (kajsa.mollersen@uit.no) 15th October 2018

% This function selects the "best" column of Z for a binary matrix
% factorization of X, measured as an increase in equal entries compared to
% mask (typically the rank k-1 approximation).

% Input:    X (n x d)    - the binary matrix that is being approximated
%           Z (n x m)    - m candidate columns (m <= n)
%           mask (n x d) - masking the elements of X not of concern

% For each column w of Z, a corresponding row vector h is constructed based
% on the similarity of w and each column of X. Then the similarity between 
% each w*h and X masked is computed. The best column w and its 
% corresponding h is returned. The column w is removed from Z. 

% Output:   w (n x 1)   - best column vector
%           h (1 x d)   - corresponding row vector
%           Z (n x m-1) - candidate columns

if ~islogical(X) || ~islogical(Z) || ~islogical(mask)
  disp('Logical, please')       % Only accept logical input
  return                        
end

[n, d] = size(X);
m = size(Z,2);                  % Number of columns in Z
H = false(m,d);                 % The row vectors that will be made

Lw = zeros(m+1,d);
for col = 1: m
  w = Z(:,col);
  for j = 1: d
    WndX = w == X(:,j);
    WndX(mask(:,j)) = false;
    Lw(col,j) = sum(WndX);
  end    
end

% Repeat for zero vector
w = false(n,1);
for j = 1: d                  % All the columns in X
  WndX = w == X(:,j);
  WndX(mask(:,j)) = false;
  Lw(m+1,j) = sum(WndX);
end

[~, idx] = max(Lw);          % w  or a vector of zeros?

idx0 = idx == m+1;    % For these columns, a zero vector is best
mask0 = ~mask;
mask0(:,idx0) = false;

for col = 1: m
  H(col,:) = idx == col;                 
end

crit = zeros(1,m);
for col = 1: m
 A = Z(:,col)*H(col,:);
 eq = A(mask0) == X(mask0);
 crit(col) = sum(eq(:));
end

[~,best_col] = max(crit);     % Which w has most 1's in the h
w = Z(:,best_col);            % That's the "best" column
sum(w)
h = H(best_col,:);            % and row

Z(:,best_col) = [];           % Remove that column from Z

