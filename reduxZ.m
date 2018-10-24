function Zredux = reduxZ(Z0, min_class, max_class, fig_nr)

Z = Z0;
n = size(Z,1);
overlap_class = n - min_class;

for t = 1: 10 % Maximum 10 iterations
  sumZcol = sum(Z,2);  
  Z(:,sumZcol > max_class) = [];
  sumZ = sum(Z,1);
  
  Zredux = [];
  used = [];
  nZ = size(Z,2);

  for k = 1 : nZ
    
    if isempty(find(used ==k,1)) % for columns that are not used yet
      
      w = Z(:, k);
  
      % These are candidates for being the same as w
      range = [sum(w)- min_class sum(w)+min_class];
      idx = find(sumZ > range(1) & sumZ < range(2));
      Lia = ismember(idx,used); 
      idx(Lia) = [];                          % Exclude the used ones
      idx(idx ==k) = [];                      % and the current column
      eq = [];
      for i = 1: length(idx) % Check if the candidates are equal
        if sum(w == Z(:,idx(i))) > overlap_class 
          eq = [eq idx(i)];
        end
      end

      if length(eq)>1
        if mod(length(eq),2) == 0 
          med_vec = median(Z(:,[eq k]),2);
        else
          med_vec = false(n,2);
          med_vec(:,1) = median(Z(:,[eq k]),2);
          med_real = median(double(Z(:,[eq k])),2);
          med_vec(med_real>0.5,2) = true;
        end
          
        if fig_nr
          figure(fig_nr), subplot(1,2,1), imagesc(Z(:,[eq k])), colormap(gray), title(k)
          subplot(1,2,2), imagesc(repmat(med_vec,1,10)), colormap(gray), title(length(eq))
          drawnow
          pause
        end
        Zredux = [Zredux med_vec];
      else
        Zredux = [Zredux w];
      end

      used = [used k eq]; 
    end
  end
  size(Zredux,2)
  if size(Zredux,2) == size(Z,2)
      t
    break
  end
  Z = logical(Zredux); % Update for the next round
end

