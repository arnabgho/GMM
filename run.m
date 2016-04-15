close all;
results = {};
fname = 'mnist';
fprintf('loading data from %s.mat\n', fname);
load(fname);


RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)))
%% do a bit of dimensionality reduction

[N D] = size(trainX);
if isfield(results, 'pca'),
  X = results.pca;
  XX = results.pcaXX;
else
  [X,evecs,evals] = PCA(trainX, 50);
  XX = X*(pinv(X)*(trainX-ones(N,1)*mean(trainX)))+ones(N,1)*mean(trainX);
  results.pca = X;
  results.pcaXX=XX;
end;
%% X is reduced, XX is reconstructed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('We have reduced dimensionality (with PCA) and will now run a simple experiment\nfor clustering into 5 clusters.  We will do this clustering 6 times\n (using different initializations) and plot\nthe Complete LL/Incomplete LL for each.\n\n');

if ~isfield(results, 'gmmcll'),
  for iter=1:6,
    [mu,pk,z,si2,CLL,ILL] = gmm(X, 5);
    results.gmmcll{iter} = CLL;
    results.gmmill{iter} = ILL;
  end;
end;
figure(1);
for iter=1:6,
  subplot(2,3,iter);
  plot(1:length(results.gmmcll{iter}),results.gmmcll{iter},'rx-', 1:length(results.gmmcll{iter}),results.gmmill{iter},'bo-');
  xlabel('iteration');
  ylabel('log likelihood');
end;
saveas(figure(1), 'fig1.jpg','jpg');
mypause


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\nNow, we will try different numbers of clusters, evaluating according\nto data likelihood and keeping the "best" results (from different initializations).\n\n');

KK=[2 5 10 15 20 25 30];
if ~isfield(results, 'clusters'),
  results.clusters = {};
  results.clustersBIC = zeros(1,length(KK));
  for ii=1:length(KK),
    fprintf('  %d clusters...', KK(ii));
    bestBIC = -Inf;
    best    = {};
    for pass=1:10,
      [mu,pk,z,si2,CLL,ILL,BIC] = gmm(X, KK(ii));
      if BIC > bestBIC,
        bestBIC = BIC;
        best.mu = mu;
        best.pk = pk;
        best.z  = z;
        best.si2 = si2;
      end;
    end;
    fprintf('\n');
  
    results.clusters{ii} = best;
    results.clustersBIC(ii) = bestBIC;
  end;
end;

figure(2);
plot(KK, results.clustersBIC, 'rx-');
xlabel('number of clusters');
ylabel('BIC');
saveas(figure(2), 'fig2.jpg','jpg');
mypause

fprintf(['\n\nFinally, we draw the means of each cluster for ' ...
         'visualization.\n\n']);
for ii=1:length(KK),
  figure(2+ii);
  % compute means
  dig = zeros(KK(ii),D);
  for k=1:KK(ii),
    dig(k,:) = sum(repmat(results.clusters{ii}.z(:,k),1,784).*trainX) / sum(results.clusters{ii}.z(:,k));
  end;
  % sort by pk
  [vv,jj] = sort(results.clusters{ii}.pk,'descend');
  draw_digits(dig(jj,:), vv);
  saveas(figure(2+ii), sprintf('fig%d.jpg',2+ii),'jpg');
end;
