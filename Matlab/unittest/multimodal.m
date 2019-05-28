K1arr = dlmread('../../unittest/multimodal/K_const/K1arr.txt',',');
K2arr = dlmread('../../unittest/multimodal/K_const/K2arr.txt',',');
logPostarr = dlmread('../../unittest/multimodal/K_const/logPostarr.txt',',');

K1arr(:,end) = [];
K2arr(:,end) = [];
logPostarr(:,end) = [];

contourf(K1arr/max(K1arr(:)), K2arr/max(K2arr(:)), logPostarr, 20)