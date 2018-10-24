function Zsets = column_set(Z)

if ~islogical(Z)
  disp('Logical, please')
  return
end

[n, nZ] = size(Z);
sumZ = sum(Z,1);
[~, idx] = sort(sumZ,'descend');
Z = Z(:,idx);
tic
Zsets = cell(1,nZ);
for k =  1: nZ
  w = Z(:,k);  
  cand = [];
  for r = 1 : 50
    w_cand = [w Z(:,cand)];
    w_or_cand = logical(sum(w_cand,2));
    pen = ones(1,nZ)*n;
    rew = zeros(1,nZ);
    crit = ones(1,nZ)*(-n);
    for i = setdiff(1:nZ,[k cand])
      pen(i) = sum(w_or_cand & Z(:,i));  
      rew(i) = sum(~logical(w_or_cand) & Z(:,i));
      crit(i) = rew(i) - pen(i);
    end
    if ~isempty(find(pen== 0,1))
      crit(pen > 0) = -n;
    end
       
    [maxcrit, idx] = max(crit);
     
    if maxcrit < 0
      Zsets{k} = [k cand];
%       if length(Zsets{k}) > 3
%         [k cand]
%       end
     % figure(55), imagesc(Z(:,[k cand])), colormap(gray)
     
      break
    else
      cand = [cand idx];
    end
    if r == 50
        pause
    end
  end
end
toc