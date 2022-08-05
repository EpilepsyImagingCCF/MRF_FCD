function new_matrix = shuffle_v2_k(old_matrix,n_fold)

[n,m] = size(old_matrix);
new_matrix = zeros(n,m);
nu_n_fold = floor(n/n_fold);
moreadd = rem(n,n_fold);
fold_matrix = [0:1:n_fold].*nu_n_fold;
fold_matrix(1,n_fold-moreadd+2:end) = fold_matrix(1,n_fold-moreadd+2:end) + [1:1:moreadd];
re = 1;

while re == 1
    re = 0;
    rand('twister', sum(10000*clock)); % random
    shuffle_number = rand(n,1);
    
    for q = 1:n
    
        [x,index] = max(shuffle_number);
        new_matrix(q,:) = old_matrix(index,:);
        shuffle_number(index) = 0;
    
    end
    
        for i = 1:n_fold  
            
              idx_test = 1+(i-1)*round(n/n_fold):i*round(n/n_fold);
              
              if max(idx_test) <= n
              if i == n_fold && max(idx_test) < n
                 idx_test = 1+(i-1)*round(n/n_fold) : n;
              end
    
              label = new_matrix(idx_test,end);
            
              if length(label) == sum(label) || sum(label) == 0
                re = 1;
              end
              end
        end
        
    
end
end


        
               
  
    
    