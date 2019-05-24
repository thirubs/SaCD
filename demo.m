%% To run this file, you need to have MATLAB tensor toolbox tool. 

%% for synthetic set parameters 
ts = 8192;%1024;%4096;%8192;%16384;%
tensor_size = [ts ts ts];
nnz = 0.0001;
tensor_rank = 100;
maxiter = 30;
X = sptenrand(tensor_size,nnz);

%% Initialization of factor matrix
N = ndims(X);
J = tensor_rank;
Uinit = cell(N,1);    
for n = 1:N 
    Uinit{n} = normalize_factor(rand(size(X,n),J),2);        
end
U=arrange(ktensor(Uinit));


%% Factorization algorithm
[X_SaCD] = SaCD(X,tensor_rank,U,maxiter); % SaCD

[X_FSaCD] = FSaCD(X,tensor_rank,U,maxiter); % FSaCD