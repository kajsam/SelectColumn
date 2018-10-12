function [w, h, Z] = select_column(X,Z,mask,K)

if ~islogical(X) || ~islogical(Z) || ~islogical(mask)
  'Logical, please'
  return
end

[n, d] = size(X);
k = size(Z,2);
cover = zeros(1,n);
H = false(k,d);
for col = 1: k
  w = Z(:,col);
  h = false(1,d);
  cov = 0;
  for j = 1 : d
    mask_col = mask(:,j);               % mask column j
    x_col = X(:,j);                     % X column j
  
    idx0 = ~mask_col;                   % idx mask = 0
    v = (x_col(idx0) & w(idx0));        % x = 1 & w = 1 
    vcov = (~x_col(idx0) & ~w(idx0));   % x = 0 & w = 0 
    idx0 = ~mask_col & ~x_col;          % idx mask = 0 & x = 0
    u = w(idx0);                        % w = 1;
    if sum(v) > sum(u)
      h(j) = 1;
      % cov = cov + sum(v) - sum(u);
      cov = cov + sum(v) + sum(vcov)-sum(u);
    else
      % cov = cov + sum(u) + sum(vcov);
    end
  end
  H(col,:) = h;
  cover(col) = cov;
end

[~,best_col] = max(cover);

w = Z(:,best_col);
h = H(best_col,:);

Z(:,best_col) = [];






    
  
  