function [Z,vecs,vals] = PCA(X,K)
  
% X is N*D input data, K is desired number of projection dimensions (assumed
% K<D).  Return values are the projected data Z, which should be N*K, vecs,
% the D*K projection matrix (the eigenvectors), and vals, which are the
% eigenvalues associated with the dimensions
  
[N D] = size(X);

if K > D,
  error('PCA: you are trying to *increase* the dimension!');
end;

% first, we have to center the data, so that the mean of each dimension
% is zero.  in other words, mean(X(:,d))==0 for all d.
mu = mean(X);
X  = X - repmat(mu,N,1);

% next, comput the covariance of the data

C = cov(X);

% compute the top K eigenvalues and eigenvectors of C... 
% hint: use 'eigs' if you're in matlab or 'eig' if you're in octave

if isMatlab,
  [vecs,vals] = eigs(C, K);
  vals = diag(vals);
else
  [vecs,vals] = eig(C);
  vals = diag(vals);
  [vals, idx] = sort(vals, 'descend');
  vecs = vecs(:,idx);
  vecs = vecs(:,1:K);
  vals = vals(1:K);
end;

Z = X*vecs;
