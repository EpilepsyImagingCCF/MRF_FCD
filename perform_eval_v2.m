function [AUC, accuracy, sensitivity, specificity, table,threshold]...
    = perform_eval_v2(label, newindex, posclass,flagCI,flagGraph)
if length(newindex) ~= length(label)
    warning('Check your inputs');
end

if flagCI == 1
    [X,Y,T,AUC] = perfcurve(label,newindex,posclass,'nboot',1000,'xvals','all');
else
    [X,Y,T,AUC] = perfcurve(label,newindex,posclass);
end

if flagGraph ==1
    plot(X,Y) 
    xlabel('1-Specificity'); ylabel('Sensitivity');
end


score = zeros(1,length(X));
for i = 1:length(X)
    score(i) = X(i,1)*X(i,1) + (1-Y(i,1))*(1-Y(i,1));
end
[c, I] = min(score);
sensitivity = Y(I,1);
specificity = 1-X(I,1);
threshold = T(I);

AC = 0;
for i=1:length(label)
    if posclass == label(i)
        AC = AC +1;
    end
end
BD = length(label)-AC;

A = round(sensitivity*AC);
D = round(specificity*BD);
C = AC-A;
B = length(label)-A-C-D;

table = [A B;C D];
accuracy = (AC*sensitivity + BD*specificity)/(A+B+C+D);
