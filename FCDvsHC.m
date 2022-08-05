%% SVM for FCD VS HC %%
n_fold = 5;
training_per = 1;
bootstrap = 100;

FCD_gaussian_long_pt = [fit_T1_pt_mean_long fit_T1_pt_std_long fit_T2_pt_mean_long fit_T2_pt_std_long];
FCD_gaussian_short_pt = [fit_T1_pt_mean_short fit_T1_pt_std_short fit_T2_pt_mean_short fit_T2_pt_std_short];

FCD_gaussian_long_hc = [fit_T1_hc_mean_long fit_T1_hc_std_long fit_T2_hc_mean_long fit_T2_hc_std_long];
FCD_gaussian_short_hc = [fit_T1_hc_mean_short fit_T1_hc_std_short fit_T2_hc_mean_short fit_T2_hc_std_short];

FCD_total_pt = [FCD_gaussian_long_pt FCD_gaussian_short_pt T1_pt_mean_data T1_pt_std_data T2_pt_mean_data T2_pt_std_data GM_pt_mean_data GM_pt_std_data WM_pt_mean_data WM_pt_std_data];
FCD_total_hc = [FCD_gaussian_long_hc FCD_gaussian_short_hc T1_hc_mean_data T1_hc_std_data T2_hc_mean_data T2_hc_std_data GM_hc_mean_data GM_hc_std_data WM_hc_mean_data WM_hc_std_data];
FCD_vs_HC = [FCD_total_pt ; FCD_total_hc];
FCD_type_label = [ones(size(FCD_total_pt,1),1);zeros(size(FCD_total_hc,1),1)];

FCD_type1 = FCD_total_pt;
FCD_type2 = FCD_total_hc;
FCD_type = FCD_vs_HC;

   
for i = 1:size(FCD_type,2)
    Fscore_num = (mean(FCD_type1(:,i)) - mean(FCD_type(:,i)))^2 + (mean(FCD_type2(:,i)) - mean(FCD_type(:,i)))^2;
    Fscore_denon = sum((FCD_type1(:,i) - mean(FCD_type1(:,i))).^2)/(size(FCD_type1,1)-1) + sum((FCD_type2(:,i) - mean(FCD_type2(:,i))).^2)/(size(FCD_type2,1)-1);
    Fscore_dd(1,i) = Fscore_num/Fscore_denon;
end
p_sub_corr = Fscore_dd;
[Fscoreo, FscoreI] = sort(p_sub_corr,'ascend');
inputdata = [FCD_type FCD_type_label];
[dx, dy] = size(inputdata);

for re = 1:bootstrap

inputdata = shuffle_v2_k(inputdata,n_fold);
    
trainingdata = inputdata(1:round(size(inputdata,1)*training_per),:);
testdata = inputdata(round(size(inputdata,1)*training_per)+1:end,:);

[dx, dy] = size(trainingdata);

ML_data = trainingdata; 
data_train_set = ML_data(:,1:dy-1); 
data_train_label = ML_data(:,dy);   
rank_SVM_index = FscoreI; 

for sub = dy-1:-1:1
cvdata = [data_train_set(:,rank_SVM_index(sub:dy-1)), data_train_label];
[cv_row, cv_col] = size(cvdata);
for i=1:n_fold
    
idx_test = 1+(i-1)*round(cv_row/n_fold):i*round(cv_row/n_fold);
    
if max(idx_test) <= cv_row
        
if i == n_fold && max(idx_test) < cv_row
idx_test = 1+(i-1)*round(cv_row/n_fold): cv_row;
end
    
cv_test_set = cvdata(idx_test,1:end-1);
cv_test_label = cvdata(idx_test,end);
tmp_data = cvdata;
tmp_data(idx_test,:) = [];
cv_training_set = tmp_data(:,1:end-1);
cv_training_label = tmp_data(:,end);
svm_structure = fitcsvm(cv_training_set,cv_training_label,'KernelFunction','linear','KernelScale','auto');
% svm_structure = fitcsvm(cv_training_set,cv_training_label,'KernelFunction','rbf','KernelScale','auto');
close all;
[result_svm_class, result_svm_regression] = predict(svm_structure,cv_test_set);  
result_svm_regression = result_svm_regression(:,1);
[~, result_acc_quad(sub,i,re), result_sen_quad(sub,i,re), result_spe_quad(sub,i,re), ~] = perform_eval_v2(cv_test_label, result_svm_regression, 0, 0, 0);    
end
end
end
tmp_ACC = result_acc_quad(:,:,re);
tmp_ACC_2(:,re) = mean(tmp_ACC,2);
tmp_SEN = result_sen_quad(:,:,re);
tmp_SEN_2(:,re) = mean(tmp_SEN,2);
tmp_SPE = result_spe_quad(:,:,re);
tmp_SPE_2(:,re) = mean(tmp_SPE,2);
end
result_mean(:,1) = mean(tmp_ACC_2,2);
result_mean(:,2) = mean(tmp_SEN_2,2);
result_mean(:,3) = mean(tmp_SPE_2,2);
result_std(:,1) = std(tmp_ACC_2');
result_std(:,2) = std(tmp_SEN_2');
result_std(:,3) = std(tmp_SPE_2');
result_mean = flipud(result_mean);
result_std = flipud(result_std);
figure(), errorbar(result_mean(:,1),result_std(:,1));
hold on
errorbar(result_mean(:,2),result_std(:,2));
errorbar(result_mean(:,3),result_std(:,3));
hold off
clearvars -except result_mean result_std 