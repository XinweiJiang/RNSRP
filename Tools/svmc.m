function [ rate,Result, optC, optG ] = svmc( x,y,xx,yy,cc, gg, tt, iscv )
%Default
%[accSvm, retYYSvm, optC, optG] = svmc(Z, Y, ZZ, YY, 0, 0, 2, 1);
addpath('./svm2/matlab');
if ~exist('iscv','var'),                  iscv = 0;  end
if ~exist('tt','var'),                      tt = 2;  end
if ~exist('gg','var'),                    gg = 1;  end
if ~exist('cc','var'),                    cc = 10;  end

if iscv
    C = [0.1,1,10,100,1000,10000]; %[1000];
    G = [0.001,0.01,0.1,1,10,100,1000]; %[0.01]
else
    C = [10];
    G = [1];
end

Acc  = [];
for cc = C
    for gg = G
        %% Training SVM
%         cmd = ['-q ' ' -c ' num2str(cc) ' -g ' num2str(gg) ' -t ' num2str(tt) ];
%         model = svmtrain(y, x, cmd);
        model = svmtrain(y, x, ['-t ' num2str(tt) ' -c ' num2str(cc) ' -g ' num2str(gg) ' -q']);
%         model = svmtrain(y, x,'-t 1 -d 2 -c 1 -g 1');

        %% Testing SVM
        [Result, accuracy, decision_values] = svmpredict(yy, xx, model);
        Acci = [size(x,1) size(xx,1) accuracy(1) cc gg];
        Acc = [Acc; Acci];
    end
end

ac = Acc(:,3);
rate = max(ac);
maxIdx = find(ac == max(ac));
if length(maxIdx) > 0
    maxIdx = maxIdx(1);
end
optC = Acc(maxIdx, 4);
optG = Acc(maxIdx, 5);

model = svmtrain(y, x, ['-t ' num2str(tt) ' -c ' num2str(optC) ' -g ' num2str(optG) ' -q']);
[Result, accuracy, decision_values] = svmpredict(yy, xx, model);

fprintf('\nBest Acc with SVM: %4.2f; optC:%6.2f; optG: %6.2f\n', rate, optC, optG);


end

