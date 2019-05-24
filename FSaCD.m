function [X_FSaCD] = FSaCD(X,J,U,maxiters)
    %% Title
    fprintf('Fast-SaCD\n');

    %% Initialization for parameters 
    N = ndims(X);
    epsilon=1e-12;
    
    %% compute Hadamard products for tensor factors
    % Hmat=sparse(ones(J,J));
    Hmat=ones(J,J);
    for n=1:N
        Hmat= Hmat .* ((U{n,1}' *U{n,1})+epsilon);
    end
    for n = 1:N
        diff{n,1} = zeros(size(U{n,1}));
    end
    
    %% iterations
    for i = 1:maxiters
        for n= 1:N
            rows_and_cols = size(U{n}); %no of rows and columns of both tensors are same
            cols = rows_and_cols(2);
            %Update Hadamard matrix for nth mode      
            Hmat = Hmat ./ ((U{n,1}' *U{n,1})+epsilon);
            %Calculate Gradient
            tmphmat = U{n,1}* Hmat;
            L{1} = norm(Hmat);
            L{2} = Hmat;
            [Unew, diff_new]=goWiter_FSaCD(X,U,L,tmphmat,cols,i,diff{n,1});         
            U{n,1}=  Unew;
            diff{n,1} = diff_new;
            %Check nonnegativity 
            U{n,1}(U{n,1}<=epsilon)=epsilon;
            %Normalization (if you need)
            if (n~=N)
              U{n,1}=normalize_factor(U{n,1},2);
            end
            %Update Hadamard matrix with updated matrix			
            Hmat = Hmat .* ((U{n,1}' *U{n,1})+epsilon); 
        end
    end
X_FSaCD = U;

