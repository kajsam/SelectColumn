function H = calculate_h(X,Z)

% Kajsa Mollersen (kajsa.mollersen@uit.no) 24th of October 2018

% Given a binary matrix X and a set of column vectors Z, this function
% constructs the corresponding set of row vectors H, such that ZH is an
% approximation to X.

% Input:        X - binary matrix
%               Z - set of columns and their combinations

if ~islogical(X) || ~islogical(Z)
  disp('Logical, please')       % Only accept logical input
  return                        
end

[n, d] = size(X);
m = size(Z,2);                  % Number of columns in Z
H = false(m,d);                 % The row vectors that will be made

Lw = zeros(m,d);
for col = 1: m
  w = Z(:,col);
  for j = 1: d
    WndX = w == X(:,j);         % counting the number of equal entries
    Lw(col,j) = sum(WndX);
  end    
end

% If a zero vector has more equal entries than any of the columns in Z,
% than the corresponding entry of h will be 0.
for j = 1: d                  
  WndX = ~X(:,j);
  Lw(m+1,j) = sum(WndX);
end

[~, idx] = max(Lw);          % Which column has the most equal entries

for col = 1: m
  H(col,:) = idx == col;                 
end