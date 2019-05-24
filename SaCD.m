function [X_out] = SaCD(X,J,U,maxiters)
    %% Title
    fprintf('SaCD\n');
	
    %% Initialization of some parameters
    N = ndims(X);
    epsilon=1e-12;
    for n = 1:N
        diff{n,1} = zeros(size(U{n,1}));
    end
		
    %% compute Hadamard products for tensor factors
    Hmat=ones(J,J);
    for n=1:N
        Hmat= Hmat .* ((U{n,1}' *U{n,1})+epsilon);
    end
    
    %% iterations
    for i = 1:maxiters
        for n= 1:N          
            factor = n;
            rows_and_cols = size(U{n}); %no of rows and columns of both tensors are same
            rows = rows_and_cols(1);
            cols = rows_and_cols(2);
            %Update Hadamard matrix for nth mode       
            Hmat = Hmat ./ ((U{n,1}' *U{n,1})+epsilon);
            tmpmat=mttkrp(X,U,n);
            mtt_time = mtt_time +toc;			
            %Calculate Gradient
            grad= -(tmpmat-(U{n,1} *Hmat));                               
            L{1} = norm(Hmat);
            L{2} = Hmat;			
			% Factor Matrix update
            [Unew, diff_new]=goWiter_SaCD(grad,L,U{n,1},cols,diff{n,1});         
            U{n,1}=  Unew;
            diff{n,1} = diff_new;
            %Check nonnegativity 
            U{n,1}(U{n,1}<=epsilon)=epsilon;		
            %Normalization (if you need)
            if (n~=N)
              U{n,1}=normalize_factor(U{n,1},1);
            end			
            %Update Hadamard matrix with updated matrix			
            Hmat = Hmat .* ((U{n,1}' *U{n,1})+epsilon); 			
        end 
    end
X_out = U;

 