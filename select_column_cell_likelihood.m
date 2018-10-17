function [w, h, Z] = select_column_cell_likelihood(X,Z,mask, alphabeta)

% Kajsa Mollersen (kajsa.mollersen@uit.no) 15th October 2018

% This function selects the "best" column of Z for a binary matrix
% factorization of X, measured as an increase in likelihood compared to
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

alpha = alphabeta(1,:);
beta = alphabeta(2,:);
t00 = alphabeta(3,:);
t01 = alphabeta(4,:);

[n, d] = size(X);
m = size(Z,2);                  % Number of columns in Z

H = false(m,d);                 % The row vectors that will be made

wask = ~mask;

Lw = zeros(m+1,d);
for col = 1: m
  w = Z(:,col);
  
  for j = 1: d                  % All the columns in X
    Lw(col,j) = -sum(t00(~w&wask(:,j))) - sum(t01(w&wask(:,j))) + sum(alpha(X(:,j)&wask(:,j)))+ sum(beta(X(:,j)&w&wask(:,j)));
  end
  
end
% Repeat for zero vector
for j = 1: d                  % All the columns in X
  Lw(m+1,j) = -sum(t00(wask(:,j))) +  sum(alpha(X(:,j)&wask(:,j)));
end
   
[~, idx] = max(Lw);          % w, ~w or a vector of zeros?

for col = 1: m
  H(col,:) = idx == col;                 % Keep h
  sumh = sum(H,2);
end
max(sumh)
% figure, imagesc(H), colormap(gray), title('H')
[~,best_col] = max(sumh);     % Which w*h has largest likelihood?
w = Z(:,best_col);              % Thats' the best column
h = H(best_col,:);              % and row

% figure, imagesc(w*h), colormap(gray), title(best_col),

Z(:,best_col) = [];             % Remove that column from Z


