function [W,diff_all] = goWiter_SaCD(GW,HH,W,k,diff_old)
lc = HH{1};
HH2 = HH{2};
s = GW/(HH2);
s = W - s;
s = max(0,s);
s = s - W;
ss = s;
diff_all = ((-1)*s.*GW-0.5*lc*s.*s);
mi=min(min(diff_all));
% ma = max(max(diff_all));
mean_def = (sum(sum(diff_all)));
mean_def_old = (sum(sum(diff_old)));		  
for p=1:k
    if (initer == 1)
        index_to_update = diff_all(:,p)>=mi;
        W(index_to_update,p) = W(index_to_update,p)+ss(index_to_update,p);
        no_of_var_updates = no_of_var_updates+sum(index_to_update);
        no_of_gradient_updates = no_of_gradient_updates+sum(index_to_update);
        % updated_variable(index_to_update,p) = updated_variable(index_to_update,p)+1; % return this if wanted to know which variables are updated based on element importance
    else
		go_to = 1;
		if mean_def > mean_def_old
			go_to = 0;
		end
		if go_to
			index_to_update =  (diff_all(:,p)<diff_old(:,p));
		else
			index_to_update =  (diff_all(:,p)>diff_old(:,p));
		end
	    W(index_to_update,p) = W(index_to_update,p)+ss(index_to_update,p);   
    end 
end


