function [Wout,diff_all] = goWiter_FSaCD(X,U,HH,tmphmat,k,initer,diff_old)
W = U{factor,1};

%% - New fast approach 
    lc = HH{1};
    HH3 = HH{2};
    N = 3;
    He = zeros(k,1);
    Z_out = cell(k,1);
    for p = 1:k
        He(p) = HH3(p,p);
        Z = cell(N,1);
        for i = [1:factor-1,factor+1:N]
            Z{i} = U{i}(:,p);
        end
        Z_out{p} = Z;
    end

    % Column-wise parallelized MTTKRP calculation and Element importance calculation
    parfor p = 1:k
        % Perform ttv multiplication
        V = double(ttv(X, Z_out{p}, -factor));
        s = -(V-tmphmat(:,p))/He(p);
        s = W(:,p)-s;
        s = max(0,s);
        s = s - W(:,p);
        ss(:,p) = s;
        diff_all(:,p) = ((-1)*s.*(-(V-tmphmat(:,p)))-0.5*lc*s.*s);
    end
    mi=min(min(diff_all));
    % ma = max(max(diff_all));\
    mean_def = (sum(sum(diff_all)));
    mean_def_old = (sum(sum(diff_old)));

    for p=1:k
        if (initer == 1)
            index_to_update = diff_all(:,p)>=mi;
            W(index_to_update,p) = W(index_to_update,p)+ss(index_to_update,p);
        else
            go_to = 1;
            if mean_def > mean_def_old
                go_to = 0;
             end
            if p == 1
                type;
            end
           if go_to
               index_to_update =  (diff_all(:,p)<diff_old(:,p));
           else
               index_to_update =  (diff_all(:,p)>diff_old(:,p));
           end
            W(index_to_update,p) = W(index_to_update,p)+ss(index_to_update,p);
        end 
    end

Wout = W;

