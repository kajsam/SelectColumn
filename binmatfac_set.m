function [Zet,H] = binmatfac_set(X, Z, Zsets, fig_nr)

% Requires:     calculate_h.m

% Given a matrix X and sets of candidate columns, this function chooses the
% best set and returns the corresponding matrix of row vectors, H. 

% Input:        X - binary matrix
%               Z - all candidate columns
%               Zsets - the entries for the sets

nZets = length(Zsets);   % The number of sets we are dealing with
FPN = zeros(1, nZets);  % The criterion for choosing the best set

for set = 1: nZets                      % For each set
  Zet = Z(:,Zsets{set});                % Extract the columns
  K = size(Zet,2);
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
  
  % figure(fig_nr+1), imagesc(Zexp), colormap(gray)
        
  H = calculate_h(X,Zexp);              % Corresponding row vectors
    
  A = logical(Zexp*H);                  % Approximation matrix
  FP = sum(sum(A & ~X));                % False positives
  FN = sum(sum(~A & X));                % False negatives
  FPN(set) = FP+FN;                     % Falses
  
  if fig_nr
    figure(fig_nr), imagesc(A), colormap(gray), title([FP FN])
   xlabel(set), ylabel(K)
  end
end  

[~,best] = min(FPN);                     % Best set

Zet = Z(:,Zsets{best});
K = size(Zet,2);
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
Hexp = calculate_h(X,Zexp);
H = Hexp(1:K,:);
for telt = K+1: size(Hexp,1)
  track = h_track{telt};
  for t = 1: length(track)
    H(track(t),:) = logical(H(track(t),:) + Hexp(telt,:));
  end
end 
  
% A = logical(Zexp*Hexp);
% FP = sum(sum(A & ~X));
% FN = sum(sum(~A & X));
% figure(fig_nr), imagesc(A), colormap(gray), title([FP FN])
% xlabel(best)

% A = logical(Zet*H);
% FP = sum(sum(A & ~X));
% FN = sum(sum(~A & X));
% figure(fig_nr+2), imagesc(A), colormap(gray), title([FP FN])
% xlabel(best)
% FPN = FP+FN;