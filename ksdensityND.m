

XS = h5read('NuclearDataFiles.hdf5',"/BoundsTMC");
en = h5read('NuclearDataFiles.hdf5',"/energies");


points = [6238, 7486];

%% Plotting marginals

[x1,y1] = ksdensity(XS(points(1),:));
[x2,y2] = ksdensity(XS(points(2),:));

figure
plot(y1,x1);
%ecdf(XS(points(1),:));
figure;
%ecdf(XS(points(2),:));
plot(y2,x2);




%%  Plotting the data marginals + the gaussians. Also show corresponding sampling

figure
scatter(XS(points(2),:),XS(points(1),:))

pdf1 = fitdist(XS(points(1),:)','kernel');
pdf2 = fitdist(XS(points(2),:)','kernel');


figure
plot(y1,pdf(pdf1,y1))
figure
plot(y2,pdf(pdf2,y2))

mean1 = mean(XS(points(1),:));
mean2 = mean(XS(points(2),:));

mu = [mean(XS(points(1),:)),mean(XS(points(2),:))];

co1 = cov([XS(points(1),:);XS(points(2),:)]');

samples = mvnrnd(mu,co1,700);


figure
scatter(XS(points(2),:),XS(points(1),:))
hold on
scatter(samples(:,2),samples(:,1))

y22 = y2(1):y2(1)/1000:y2(end)+y2(end)/100;

marg1 = normpdf(y1,mu(1),sqrt(co1(1,1)));
marg2 = normpdf(y22,mu(2),sqrt(co1(2,2)));

figure
plot(y1,pdf(pdf1,y1))
hold on
plot(y1,marg1);

figure
plot(y2,pdf(pdf2,y2))
hold on
plot(y22,marg2);


%% Plotting marginals plus data comparison using 


%[icdf1,y1] = ksdensity(XS(points(1),:),'function','icdf');
%[icdf2,y2] = ksdensity(XS(points(2),:),'function','icdf');

pdf1 = fitdist(XS(points(1),:)','kernel');
pdf2 = fitdist(XS(points(2),:)','kernel');

corr1 =  corr([XS(points(1),:);XS(points(2),:)]');
corr1 = corr1(1,2);

%corr = 0.8;             %Select a Pearson correlation

a1 = rand(1,700);       %Independant standard uniform
a2 = rand(1,700);

w1 = norminv(a1);   %Transform to standard normal space
w2 = norminv(a2);   %Inverse cdf method

z1 = w1;
z2 = corr1*w1+sqrt(1-corr1^2)*w2; %Induce correlation in standard normal space using 

%figure;
%plot(z1,z2,'*');

%figure;
%plot(normcdf(z1),normcdf(z2),'*');

samples1 = icdf(pdf1,normcdf(z1));
samples2 = icdf(pdf2,normcdf(z2));

figure
scatter(XS(points(2),:),XS(points(1),:))
hold on
scatter(samples2,samples1)


%figure
%plot(z1,z2 ,'*');
%figure
%plot(samples1,samples2,'*');




%% Plotting marginals plus data comparison assuming independance

% 
% %[icdf1,y1] = ksdensity(XS(points(1),:),'function','icdf');
% %[icdf2,y2] = ksdensity(XS(points(2),:),'function','icdf');
% 
% pdf1 = fitdist(XS(points(1),:)','kernel');
% pdf2 = fitdist(XS(points(2),:)','kernel');
% 
% %corr1 =  corr([XS(points(1),:);XS(points(2),:)]');
% corr1 = 0;
% 
% %corr = 0.8;             %Select a Pearson correlation
% 
% a1 = rand(1,700);       %Independant standard uniform
% a2 = rand(1,700);
% 
% w1 = norminv(a1);   %Transform to standard normal space
% w2 = norminv(a2);   %Inverse cdf method
% 
% z1 = w1;
% z2 = corr1*w1+sqrt(1-corr1^2)*w2; %Induce correlation in standard normal space using 
% 
% %figure;
% %plot(z1,z2,'*');
% 
% %figure;
% %plot(normcdf(z1),normcdf(z2),'*');
% 
% samples1 = icdf(pdf1,normcdf(z1));
% samples2 = icdf(pdf2,normcdf(z2));
% 
% figure
% scatter(XS(points(2),:),XS(points(1),:))
% hold on
% scatter(samples2,samples1)


%figure
%plot(z1,z2 ,'*');
%figure
%plot(samples1,samples2,'*');



%% Trying to bootsrap the data



pdf1 = fitdist(XS(points(1),:)','kernel');
pdf2 = fitdist(XS(points(2),:)','kernel');

pseudoSamples1 = cdf(pdf1,XS(points(1),:));
pseudoSamples2 = cdf(pdf2,XS(points(2),:));

a = empcopulapdf([pseudoSamples1;pseudoSamples2]',0.001,100,'betak-matlab');

outSamples = empcopularnd(a,5000);


figure
scatter(pseudoSamples1,pseudoSamples2)
hold on
scatter(outSamples(:,1),outSamples(:,2))

Bootsraped1 = icdf(pdf1, outSamples(:,1));
Bootsraped2 = icdf(pdf2, outSamples(:,2));


figure
scatter(XS(points(2),:),XS(points(1),:))
hold on
scatter(Bootsraped2,Bootsraped1)


cop1 = empcopulacdf([pseudoSamples1;pseudoSamples2]',100,'deheuvels');
[X,Y] = meshgrid([0.01:0.01:1],[0.01:0.01:1]);
figure
mesh(X,Y,cop1)
figure
mesh(X,Y,a)


% 
% %[icdf1,y1] = ksdensity(XS(points(1),:),'function','icdf');
% %[icdf2,y2] = ksdensity(XS(points(2),:),'function','icdf');
% 
% pdf1 = fitdist(XS(points(1),:)','kernel');
% pdf2 = fitdist(XS(points(2),:)','kernel');
% 
% %corr1 =  corr([XS(points(1),:);XS(points(2),:)]');
% corr1 = 0;
% 
% %corr = 0.8;             %Select a Pearson correlation
% 
% a1 = rand(1,700);       %Independant standard uniform
% a2 = rand(1,700);
% 
% w1 = norminv(a1);   %Transform to standard normal space
% w2 = norminv(a2);   %Inverse cdf method
% 
% z1 = w1;
% z2 = corr1*w1+sqrt(1-corr1^2)*w2; %Induce correlation in standard normal space using 
% 
% %figure;
% %plot(z1,z2,'*');
% 
% %figure;
% %plot(normcdf(z1),normcdf(z2),'*');
% 
% samples1 = icdf(pdf1,normcdf(z1));
% samples2 = icdf(pdf2,normcdf(z2));
% 
% figure
% scatter(XS(points(2),:),XS(points(1),:))
% hold on
% scatter(samples2,samples1)
% 
% 
% %figure
% %plot(z1,z2 ,'*');
% %figure
% %plot(samples1,samples2,'*');










