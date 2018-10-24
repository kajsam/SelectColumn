function [w, h, Z] = select_column_set(X,Z,mask)

% Kajsa Mollersen (kajsa.mollersen@uit.no) 23rd of October 2018

if ~islogical(X) || ~islogical(Z) || ~islogical(mask)
  disp('Logical, please')       % Only accept logical input
  return                        
end

[n, d] = size(X);
m = size(Z,2);                  % Number of columns in Z
% H = false(m,d);                 % The row vectors that will be made

wask = ~mask;

Lw = zeros(m,d);
for col = 1: m
  w = Z(:,col);
  for j = 1: d
    mask_col = wask(:,j);
    WndX = w(mask_col) == X(mask_col,j); 
    Lw(col,j) = sum(WndX);
  end    
end

% Repeat for zero vector
for j = 1: d                  % All the columns in X
  WndX = ~X(:,j);
  WndX(mask(:,j)) = false;
  Lw(m+1,j) = sum(WndX);
end

[~, idx] = max(Lw);          % w  or a vector of zeros?

% idx0 = idx == m+1;    % For these columns, a zero vector is best
% mask0 = ~mask;
% mask0(:,idx0) = false;
% 
% for col = 1: m
%   H(col,:) = idx == col;                 
% end
% 
% crit = zeros(1,m);
% for col = 1: m
%  A = Z(:,col)*H(col,:);
%  eq = A(mask0) == X(mask0);
%  crit(col) = sum(eq(:));
% end
% 
% [~,best_col] = max(crit);     % Which wh is most similar to X
w = Z(:,1);            % That's the "best" column
h = idx == 1;       

Z(:,1) = [];

