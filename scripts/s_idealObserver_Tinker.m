
n = 12;
numtrials = 200;
dozscore = true; 

f1  = 2;
f2  = f1;
ph1 = 0;
ph2 = pi/2;

mn = 100;
c  = 5;
template =  mn+c*[sin(f1*(1:n)/n*2*pi+ph1); sin(f2*(1:n)/n*2*pi+ph2)];

classes = zeros(numtrials, 2, 'logical');
classes(1:numtrials/2,1)       = 1;
classes(numtrials/2+1:end, 2)  = 1;

pred = classes * template ;


data = poissrnd(pred);

figure(1), subplot(1,2,1); imagesc(pred, [0 2*mn]), title('Templates')
subplot(1,2,2); imagesc(data, [0 2*mn]), title('Data')

%
classifierOutput = zeros(numtrials,3);
classifierOutput(:,1) = data * diff(template)';
classifierOutput(:,2) = data * diff(log(template))';
classifierOutput(:,3) = diff(corr(data', template'),[],2);


cv = fitcsvm(data, classes*[2;1], 'Standardize', true, 'KernelFunction', 'linear', 'kFold', 10);
classLoss = kfoldLoss(cv);
P = (1-classLoss) * 100;

if dozscore, classifierOutput = zscore(classifierOutput); end

classification = classifierOutput>0;

classificationSVM   = kfoldPredict(cv);

figure(2),clf 
for ii = 1:3
    subplot(3,1,ii)
    if dozscore
        h1=histogram(classifierOutput(classes(:,1),ii), -3:.1:3); hold on;
        h2=histogram(classifierOutput(classes(:,2),ii), -3:.1:3);
    else
        h1=histogram(classifierOutput(classes(:,1),ii)); hold on;
        h2=histogram(classifierOutput(classes(:,2),ii));
        h2.BinWidth = h1.BinWidth;
    end
    
end
 
figure(3), clf
for ii = 1:3
    for jj = 1:3
        subplot(3,3,(ii-1)*3+jj)
        scatter(classifierOutput(:,ii), classifierOutput(:,jj), '.');
        axis([-3 3 -3 3]); axis square;
    end        
end


figure(4); clf
for ii = 1:3
    gscatter(1:numtrials, classification(:,ii)+ii*2, classes*[2;1])
    hold on;
end
gscatter(1:numtrials, classificationSVM+8, classes*[2;1])

