function Zredux = reduxZoverlap(Z, min_class)

if ~islogical(Z)
  disp('Logical, please')
  return
end

[n, nZ] = size(Z);
sumZ = sum(Z,1);
[~, idx] = sort(sumZ,'descend');
Z = Z(:,idx);
% figure, subplot(1,3,1), imagesc(Z), colormap(gray)
keep = false(1,nZ);
for k = 1 : nZ % Check the columns one by one
  cand = true(1,nZ);
  % If the other column overlaps the 0's in the column in question, it is
  % not a representative
  for i = k : nZ
    if sum(Z(~Z(:,k),i)) > min_class
      cand(i) = false;
    end
  end
  cand(1:k) = false;
  
  Zcand = Z(:,cand);
  repr = logical(sum(Zcand,2));
%   subplot(1,3,2), imagesc([Z(:,k) repr (Z(:,k) & repr)]), colormap(gray), 
%  title('Original Representation Sum')

  if sum(Z(:,k)) - sum(Z(:,k) & repr) >  min_class
    keep(k) = true;
  end
end

Zredux = Z(:,keep);

%  subplot(1,3,3), imagesc(Zredux), colormap(gray)
 

