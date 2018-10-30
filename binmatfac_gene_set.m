function [Zet,H] = binmatfac_gene_set(X, Z, Zsets, pi, fig_nr)

% Requires:     calculate_h.m

% Given a matrix X and sets of candidate columns, this function chooses the
% best set and returns the corresponding matrix of row vectors, H. 

% Input:        X - binary matrix
%               Z - all candidate columns
%               Zsets - the entries for the sets

[n,d] = size(X)
logpi0 = log(pi(1,:));
logpi0compl = log(ones(d,1)- pi(1,:));
logpi1 = log(pi(2,:));
logpi1compl = log(ones(d,1)- pi(2,:));

nZets = length(Zsets)   % The number of sets we are dealing with
FPN = zeros(1, nZets);  % The criterion for choosing the best set
[n,d] = size(X);

for set = 1: nZets                      % For each set
  Zet = Z(:,Zsets{set});                % Extract the columns
  K = size(Zet,2)
  
  figure(fig_nr),subplot(1,2,1), imagesc(Zet), colormap(gray)
  if K > 7 % Do classical greedy row construction
    miss = n*d;
    W = false(n,K); 
    H = false(K,d);
    mask = false(size(X));
    tic
    for  it = 1: 10
      Zet = Z(:,Zsets{set});
      for k = 1: K
        [w, h, Zet] = select_column_set(X,Zet,mask);
        W(:,k) = w;
        H(k,:) = h;
        mask = mask | w*h;
      end  

      FP = sum(sum(mask & ~X));
      FN = sum(sum(~mask & X));
   
      if miss <= FP+FN
        figure(fig_nr), imagesc(mask), colormap(gray), title([FP FN])
        xlabel(set),  drawnow
        break
      else
        miss = FP+FN;
        FPN(set) = miss;
      end
    end
  else
    
    % figure(fig_nr+2), imagesc(Zet), colormap(gray)
    Zexp = Zet;
    % Expand the set to all combinations
    for k = 2 : K-1
      c = combnk(1:K,k);
      for i = 1: size(c,1)
        Zexp = [Zexp logical(sum(Zet(:,c(i,:)),2))];
      end
    end
    Zexp = [Zexp logical(sum(Zet(:,1:K),2))];
  
    H = calculate_h_gene(X,Zexp,pi);              % Corresponding row vectors
   
    A = logical(Zexp*H);                  % Approximation matrix
    
    loglik = zeros(1,d);
    for j = 1: d
      loglik(j) = sum(A(:,j) == X(:,j))*logpi1(j) + sum(~A(:,j) == X(:,j))*logpi0(j) ...
        + sum(A(:,j) == ~X(:,j))*logpi1compl(j) + sum(~A(:,j) == ~X(:,j))*logpi0compl(j);
    end  
    
    FPN(set) = sum(loglik);                     % Falses
  
    if fig_nr
      figure(fig_nr), subplot(1,2,2), imagesc(A), colormap(gray), title(FPN(set))
      xlabel(set), ylabel(K)
    end
  end
end

[~,best] = max(FPN)                     % Best set

Zet = Z(:,Zsets{best});
K = size(Zet,2);

if K > 7 % Do classical greedy row construction
  miss = n*d;
  W = false(n,K); 
  H = false(K,d);
  mask = false(size(X));
  tic
  for  it = 1: 10
    Zet = Z(:,Zsets{set});
    for k = 1: K
      [w, h, Zet] = select_column_set(X,Zet,mask);
      W(:,k) = w;
      H(k,:) = h;
      mask = mask | w*h;
    end  

    FP = sum(sum(mask & ~X));
    FN = sum(sum(~mask & X));
   
    if miss <= FP+FN
      figure(fig_nr), imagesc(mask), colormap(gray), title([FP FN])
      xlabel(set),  drawnow
      break
    else
      miss = FP+FN;
      FPN(set) = miss;
    end
  end
else
  Zexp = Zet;
  h_track = cell(1);
  for k = 1:K
    h_track{k} = K;
  end
  telt = K;
  % Expand the set to all combinations
  for k = 2 : K-1
    c = combnk(1:K,k);
    for i = 1: size(c,1)
      telt = telt+1;
      Zexp = [Zexp logical(sum(Zet(:,c(i,:)),2))];
      h_track{telt} = c(i,:);
    end
  end
  Zexp = [Zexp logical(sum(Zet(:,1:K),2))];
  h_track{telt+1} = 1:K;
      
  % figure(fig_nr+1), imagesc(Zet), colormap(gray), 
  % title('Candidate columns'), xlabel(best)
  Hexp = calculate_h_gene(X,Zexp,pi);
  H = Hexp(1:K,:);
  for telt = K+1: size(Hexp,1)
    track = h_track{telt};
    for t = 1: length(track)
      H(track(t),:) = logical(H(track(t),:) + Hexp(telt,:));
    end
  end 
end
  