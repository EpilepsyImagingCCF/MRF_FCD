%% SVM for FCD subtyping %%
n_fold = 5;
training_per = 1;
bootstrap = 100;

FCD_gaussian_long_pt = [fit_T1_pt_mean_long fit_T1_pt_std_long fit_T2_pt_mean_long fit_T2_pt_std_long];
FCD_gaussian_short_pt = [fit_T1_pt_mean_short fit_T1_pt_std_short fit_T2_pt_mean_short fit_T2_pt_std_short];
FCD_gaussian_long_hc = [fit_T1_hc_mean_long fit_T1_hc_std_long fit_T2_hc_mean_long fit_T2_hc_std_long];
FCD_gaussian_short_hc = [fit_T1_hc_mean_short fit_T1_hc_std_short fit_T2_hc_mean_short fit_T2_hc_std_short];

fcd1=[2,3,4,7,8,15,16,17,19,20,22]; fcd2=[1,5,6,9,10,11,12,13,14,18,21,23,25,24,26,27,28,29];
fcd1_ab = [FCD_gaussian_long_pt(fcd1,:) FCD_gaussian_short_pt(fcd1,:) T1_pt_mean_data(fcd1,:) T1_pt_std_data(fcd1,:) T2_pt_mean_data(fcd1,:) T2_pt_std_data(fcd1,:) GM_pt_mean_data(fcd1,:) GM_pt_std_data(fcd1,:) WM_pt_mean_data(fcd1,:) WM_pt_std_data(fcd1,:)];
fcd2_ab = [FCD_gaussian_long_pt(fcd2,:) FCD_gaussian_short_pt(fcd2,:) T1_pt_mean_data(fcd2,:) T1_pt_std_data(fcd2,:) T2_pt_mean_data(fcd2,:) T2_pt_std_data(fcd2,:) GM_pt_mean_data(fcd2,:) GM_pt_std_data(fcd2,:) WM_pt_mean_data(fcd2,:) WM_pt_std_data(fcd2,:)];
fcd1_hc_ab = [FCD_gaussian_long_hc(fcd1,:) FCD_gaussian_short_hc(fcd1,:) T1_hc_mean_data(fcd1,:) T1_hc_std_data(fcd1,:) T2_hc_mean_data(fcd1,:) T2_hc_std_data(fcd1,:) GM_hc_mean_data(fcd1,:) GM_hc_std_data(fcd1,:) WM_hc_mean_data(fcd1,:) WM_hc_std_data(fcd1,:)];
fcd2_hc_ab = [FCD_gaussian_long_hc(fcd2,:) FCD_gaussian_short_hc(fcd2,:) T1_hc_mean_data(fcd2,:) T1_hc_std_data(fcd2,:) T2_hc_mean_data(fcd2,:) T2_hc_std_data(fcd2,:) GM_hc_mean_data(fcd2,:) GM_hc_std_data(fcd2,:) WM_hc_mean_data(fcd2,:) WM_hc_std_data(fcd2,:)];
FCD_type1 = fcd1_ab./fcd1_hc_ab*100;
FCD_type2 = fcd2_ab./fcd2_hc_ab*100;

% fcd2a=[9,10,11,13,21,23,24,28]; fcd2b=[1,5,6,12,14,18,25,26,27,29];
% fcd2a_ab = [FCD_gaussian_long_pt(fcd2a,:) FCD_gaussian_short_pt(fcd2a,:) T1_pt_mean_data(fcd2a,:) T1_pt_std_data(fcd2a,:) T2_pt_mean_data(fcd2a,:) T2_pt_std_data(fcd2a,:) GM_pt_mean_data(fcd2a,:) GM_pt_std_data(fcd2a,:) WM_pt_mean_data(fcd2a,:) WM_pt_std_data(fcd2a,:)];
% fcd2b_ab = [FCD_gaussian_long_pt(fcd2b,:) FCD_gaussian_short_pt(fcd2b,:) T1_pt_mean_data(fcd2b,:) T1_pt_std_data(fcd2b,:) T2_pt_mean_data(fcd2b,:) T2_pt_std_data(fcd2b,:) GM_pt_mean_data(fcd2b,:) GM_pt_std_data(fcd2b,:) WM_pt_mean_data(fcd2b,:) WM_pt_std_data(fcd2b,:)];
% fcd2a_hc_ab = [FCD_gaussian_long_hc(fcd2a,:) FCD_gaussian_short_hc(fcd2a,:) T1_hc_mean_data(fcd2a,:) T1_hc_std_data(fcd2a,:) T2_hc_mean_data(fcd2a,:) T2_hc_std_data(fcd2a,:) GM_hc_mean_data(fcd2a,:) GM_hc_std_data(fcd2a,:) WM_hc_mean_data(fcd2a,:) WM_hc_std_data(fcd2a,:)];
% fcd2b_hc_ab = [FCD_gaussian_long_hc(fcd2b,:) FCD_gaussian_short_hc(fcd2b,:) T1_hc_mean_data(fcd2b,:) T1_hc_std_data(fcd2b,:) T2_hc_mean_data(fcd2b,:) T2_hc_std_data(fcd2b,:) GM_hc_mean_data(fcd2b,:) GM_hc_std_data(fcd2b,:) WM_hc_mean_data(fcd2b,:) WM_hc_std_data(fcd2b,:)];
% FCD_type1 = fcd2a_ab./fcd2a_hc_ab*100;
% FCD_type2 = fcd2b_ab./fcd2b_hc_ab*100;

FCD_type1_label = zeros(size(FCD_type1,1),1);
FCD_type2_label = ones(size(FCD_type2,1),1);
FCD_type = [FCD_type1 ; FCD_type2];
FCD_type_label = [FCD_type1_label ; FCD_type2_label];

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
svm_structure = fitcsvm(cv_training_set,cv_training_label,'KernelFunction','rbf','KernelScale','auto');
%     svm_structure = fitcsvm(cv_training_set,cv_training_label,'KernelFunction','linear','KernelScale','auto');
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