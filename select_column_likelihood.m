function [w, h, Z] = select_column_likelihood(X,Z,mask)

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

[n, d] = size(X);
m = size(Z,2);                  % Number of columns in Z
H = false(m,d);                 % The row vectors that will be made
loglik = zeros(1,m);            % The measure of increase in likelihood

% H_compl = false(size(H));     % Comparing to the complement
% Z_compl = ~Z;
% loglik_compl = zeros(1,K);

for col = 1: m
  w = Z(:,col);
 
  % Compare w to all columns in X
  W = repmat(w,1,d);            % For each column, compared to all columns in X
  WndX = W == X;                % Which entries are equal
  WndX(mask) = false;           % These will not count
  Lw = sum(WndX,1);             % Sum up column-wise
  
  %   % Repeat for the complement          
  %   W = repmat(~w,1,d);         
  %   WndX = W == X;   
  %   WndX(mask) = false;
  %   Lw_c = sum(WndX,1);    
  
  % Repeat for zero vector
  W = false(n,d);           
  WndX = W == X;   
  WndX(mask) = false;
  L0 = sum(WndX,1);  
    
  L = [Lw; L0]; %Lw_c;          % Which has more similar entries:
  [~, idx] = max(L);            % w or a vector of zeros?
  h = idx==1;                   % If it's w, the entry of h is 1, else 0
  H(col,:) = h;                 % Keep h
  
  eq = X == w*h;                        % Equal entries for W and w*h
  loglik(col) = sum(sum(eq(~mask)));    % Masked
  
  %   h_compl = idx==2;
  %   eq = X == ~w*h_compl;
  %   loglik_compl(col) = sum(sum(eq(~mask)));
  %   H_compl(col,:) = h_compl;
      
end

[~,best_col] = max(loglik);     % Which w*h has most equal entries?
w = Z(:,best_col);              % Thats' the best column
h = H(best_col,:);              % and row

Z(:,best_col) = [];             % Remove that column from Z

% [~,best_compl] = max(loglik_compl);
% 
% w_compl = Z_compl(:,best_compl);
% h_compl = H_compl(best_compl,:);
% 
% if max(loglik_compl)> max(loglik)
%     
%    'hey'
% end


