function [w, h, Z] = select_column_likelihood(X,Z,mask)

if ~islogical(X) || ~islogical(Z) || ~islogical(mask)
  'Logical, please'
  return
end

[n, d] = size(X);
K = size(Z,2)
H = false(K,d);
% H_compl = false(size(H));
% Z_compl = ~Z;

loglik = zeros(1,K);
% loglik_compl = zeros(1,K);
for col = 1: K
  w = Z(:,col);
 
  % Compare w to all columns in X
  % The conditional probability

  W = repmat(w,1,d);              % For each column, compared to all columns in X
  WndX = W == X;                   % Which entries are equal
  WndX(mask) = false;   % These will not count
  Lw = sum(WndX,1);    % Sum up column-wise
  
%   % Repeat for the complement          
%   W = repmat(~w,1,d);         
%   WndX = W == X;   
%   WndX(mask) = false;
%   Lw_c = sum(WndX,1);    
  
  % Reapeat for zero vector
  W = false(n,d);           
  WndX = W == X;   
  WndX(mask) = false;
  L0 = sum(WndX,1);  
    
  L = [Lw; L0]; %Lw_c; 
    
  [~, idx] = max(L);
  
  h = idx==1;  
  eq = X == w*h;
  loglik(col) = sum(sum(eq(~mask)));
  H(col,:) = h;
  
%   h_compl = idx==2;
%   eq = X == ~w*h_compl;
%   loglik_compl(col) = sum(sum(eq(~mask)));
%   H_compl(col,:) = h_compl;
  
end
[~,best_col] = max(loglik);

w = Z(:,best_col);
h = H(best_col,:);

% [~,best_compl] = max(loglik_compl);
% 
% w_compl = Z_compl(:,best_compl);
% h_compl = H_compl(best_compl,:);
% 
% if max(loglik_compl)> max(loglik)
%     
%    'hey'
% end

Z(:,best_col) = [];

