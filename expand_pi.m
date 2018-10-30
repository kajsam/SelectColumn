function [Wexp, Hexp, notWexp, notHexp] = expand_pi(W,X,pi, fig_nr)

% Expand the set to all combinations
[n,d] = size(X);
K = size(W,2)
Zexp = W;
for k = 2 : K-1
  c = combnk(1:K,k);
  for i = 1: size(c,1)
    Zexp = [Zexp logical(sum(W(:,c(i,:)),2))];
  end
end
Wexp = [Zexp logical(sum(W(:,1:K),2))];
Hexp = calculate_h_gene(X,Wexp,pi);

figure(fig_nr+1), subplot(1,3,1), imagesc(logical(Wexp*Hexp)), colormap(gray)  
title(strcat('Initial A, K = ', num2str(size(W,2))))
subplot(1,3,2), imagesc(Wexp), title('Initial columns')
subplot(1,3,3), imagesc(Hexp), title('Initial rows')

% Expand the set to all combinations
Zexp = ~W;
for k = 2 : K-1
  c = combnk(1:K,k);
  for i = 1: size(c,1)
    Zexp = [Zexp ~logical(sum(W(:,c(i,:)),2))];
  end
end
notWexp = [~Wexp ones(n,1)]; 
notHexp = [Hexp; ~sum(Hexp)];

figure(fig_nr+2)
subplot(1,3,2), imagesc(notWexp), title('Initial columns')
subplot(1,3,3), imagesc(notHexp), title('Initial rows')
subplot(1,3,1), imagesc(notWexp*notHexp), title('Initial rows')