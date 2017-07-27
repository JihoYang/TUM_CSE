height = [47;93;53;45;67;42];
weight = [15;35;15;10;27;10];

A = [height weight];

[U,S,V] = svd(A,0); sigma = diag(S);

size = sigma(1)*U(:,1);

pca_height = size*V(1,1);
pca_weight = size*V(2,1);

figure(1)
bar([height,pca_height]); title('height');
legend('data','pca')

figure(2)
bar([weight,pca_weight]); title('weight');
legend('data','pca')