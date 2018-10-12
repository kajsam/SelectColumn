function [w, h, Z] = select_column(X,Z,mask)

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
    mask_col = mask(:,j);
    x_col = X(:,j);
  
    idx0 = ~mask_col;
    v = (x_col(idx0) & w(idx0));
    idx0 = ~mask_col & ~x_col;
    u = w(idx0);
    if sum(v) > sum(u)
      h(j) = 1;
      cov = cov + sum(v)-sum(u);
    end
  end
  H(col,:) = h;
  cover(col) = cov;
end

[~,best_col] = max(cover);

w = Z(:,best_col);
h = H(best_col,:);

Z(:,best_col) = [];






    
  
  