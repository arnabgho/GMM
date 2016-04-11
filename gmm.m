function [mu,pk,z,si2,CLL,ILL,BIC] = gmm(X,K)
% input X is N*D data, K is number of clusters desired; numIter is the
% maximum number of EM iterations to run.  the output is: mu, which is K*D, the
% coordinates of the means; pk, which is K*1 and represents the cluster
% proportions; z, which is N*K, stores in positions z(n,k) the probability
% that the nth data point belongs to cluster k, specifying the cluster
% associated with each data point; si2 is the estimated (shared) variance of
% the data; CLL is a plot of the complete log likelihoods and ILL is a plot
% of the incomplete log likelihoods. BIC is the Bayesian Information
% Criterion (smaller BIC is better)
  
[N D] = size(X);

if K >= N,
  error('you are trying to make too many clusters!');
end;

numIter = 20;  % maximum number of iterations to run

% initialize the log likelihoods
CLL = zeros(numIter,1);
ILL = zeros(numIter,1);

% initialize si2 dumbly
si2 = 1;

% initialize pk uniformly
pk = ones(K,1) / K;

% we initialize the means totally randomly
mu = randn(K,D);
z = zeros(N,K);
for iter=1:numIter,
  % in the first step, we do assignments: each point is probabilistically
  % assigned to each center
  
  for n=1:N,
    den=0;
    for k=1:K,
      % compute z(n,k) = log probability that the nth data point belongs
      % to cluster k
      % z(n,k) = ((pk(k,:))*normpdf(X(n,:),mu(k,:),sqrt(si2)))
      % disp(size(X(n,:)))
      % disp("------------------------------\n")
      % disp(size(mu(k,:)))
      % disp("-------------------------------\n")
      z(n,k)=mvnpdf(X(n,:),mu(k,:),si2*eye(D));
      z(n,k)=z(n,k)*pk(k,:);
      den=den+z(n,k);
    end;
    z(n,:)=z(n,:)/den;
    % turn log probabilities into actual probabilities
    maxZ   = max(z(n,:));
    z(n,:) = exp(z(n,:) - maxZ - log(sum(exp(z(n,:) - maxZ))));
  end;
  
  Num_eff=zeros(K,1);

  for k=1:K
    for n=1:N
      Num_eff(k,1)=Num_eff(k,1)+z(n,k);
    end;
  end;    
  % re-estimate pk
  pk = Num_eff/N;
  
  % re-estimate the means
  % mu = ???
  mu=zeros(K,D);
  for k=1:K
    for n=1:N
      mu(k,:)=mu(k,:)+z(n,k)*X(n,:);
    end;
    mu(k,:)=mu(k,:)/Num_eff(k,1);
  end;    
  
  % re-estimate the variance
  % si2 = ???
  % Computing si2
  sum_sq=0;
  for n=1:N
    for k=1:K
      sum_sq=sum_sq+z(n,k)*(norm(X(n,:)-mu(k,:))^2);
    end; 
  end;   
  si2=sum_sq/(N*D);    
  
% compute *complete* and *incomplete* log likelihoods
cll = 0;
ill = 0;
for n=1:N
    ill_tmp = zeros(1, K);
    for k=1:K
        ill_tmp(k) = ill_tmp(k) + log(pk(k)) -(0.5/si2)*(X(n,:)-mu(k,:))*(X(n,:)-mu(k,:))' - 0.5*D*log(2*pi) - 0.5*D*log(si2);
        cll = cll + z(n,k)*ill_tmp(k);
    end
    % internal sum of log likelihood (adding probabilites)
    ill = ill + logsumexp(ill_tmp,2);
    % external sum of log likelihood (adding logarithms)
end
  % store values
  CLL(iter) = cll;
  ILL(iter) = ill;
end;

BIC = -2*cll+ (numel(pk)+numel(mu)+1)*log(N);
