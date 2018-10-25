function Zsets = column_set(Z, delete_noise)

% Creates a set of columns from Z for each single column of Z

% Input:    Z : candidate columns

% Output:   Zsets : indexes for each set

if ~islogical(Z)
  disp('Logical, please')
  return
end

[n, nZ] = size(Z);
Zsets = cell(1,nZ);
for k =  1: nZ
    
  w = Z(:,k);  
  cand = [];
  for r = 1 : 100
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
     
if delete_noise
     del = [];
     for i = 1: length(Zsets{k})
       if sum(Z(:,Zsets{k}(i)) & logical(sum(Z(:,setdiff(Zsets{k},Zsets{k}(i))),2))) == sum(Z(:,Zsets{k}(i)))
         del = [del i];
       end
       
     end
%      if ~isempty(del)
     %   figure(54), imagesc(Z(:,Zsets{k})), colormap(gray)
       Zsets{k}(del) = [];
 %     figure(53), imagesc(Z(:,Zsets{k})), colormap(gray), xlabel(k)
%   end
end
      break
    else
      cand = [cand idx];
    end
    if r == 100
        pause
    end
  end
end