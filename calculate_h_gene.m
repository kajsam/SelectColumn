function H = calculate_h_gene(X,Z,pi)

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
logpi0 = log(pi(1,:));
logpi0compl = log(ones(d,1)- pi(1,:));
logpi1 = log(pi(2,:));
logpi1compl = log(ones(d,1)- pi(2,:));
for col = 1: m
  w = Z(:,col);
  for j = 1: d
    Lw(col,j) = sum(w == X(:,j))*logpi1(j) + sum(~w == X(:,j))*logpi0(j) ...
        + sum(w == ~X(:,j))*logpi1compl(j) + sum(~w == ~X(:,j))*logpi0compl(j);
  end    
end

% If a zero vector has more equal entries than any of the columns in Z,
% than the corresponding entry of h will be 0.
w = false(size(w));
for j = 1: d
  Lw(m+1,j) = sum(w == X(:,j))*logpi1(j) + sum(~w == X(:,j))*logpi0(j) ...
      + sum(w == ~X(:,j))*logpi1compl(j) + sum(~w == ~X(:,j))*logpi0compl(j);
end    

[~, idx] = max(Lw);          % Which column has the most equal entries

for col = 1: m
  H(col,:) = idx == col;                 
end