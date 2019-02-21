clear all; close all;

a = [0:0.01:1];
b = [0:0.01:1];
[X,Y]=meshgrid(a,b);

%% Freche bounds + independant

op= @(x,y) max(x+y-1,0);
indep= @(x,y) x.*y;
perf= @(x,y) min(x,y);

ZZ = indep(X,Y);
figure;
meshc(X,Y,ZZ)

%% Gaussian Copula

mu = [0,0];
sigma=[1,0.5;0.5,1];

Gau=@(x,y)mvncdf([norminv(x,mu(1),sqrt(sigma(1))),norminv(y,mu(2),sqrt(sigma(4)))],mu,sigma);


ZZZ = zeros(size(X));

for i=1:length(X)
    for j=1:length(Y)
        ZZZ(i,j) = Gau(X(i,j),Y(i,j));
    end
end

figure
meshc(X,Y,ZZZ);


%% Frank copula

% s>0; 0 for perfect, 1 for indep, inf for oposite
s =200;
C = @(x,y) log(1+(s.^x-1).*(s.^y-1)./(s-1))/log(s);
Z = C(X,Y);
%figure;
%meshc(X,Y,Z);

%% Clayton Copula

% t>-1; -1 for opisite, 0 for indep and inf for perfect
t = -1;
Cla = @(x,y) max((x.^(-t)+y.^(-t)-1).^(-1/t),0);
Z2 = Cla(X,Y);
%figure
%meshc(X,Y,Z2);



%% Generating joint distribution using marginals a copula

a = [-5:0.01:5];
b = [-5:0.01:5];

[X,Y]=meshgrid(a,b);

Marginal1 = @(x) normcdf(x,0,1);
Marginal2 = @(x) normcdf(x,0,1);

Joint = @ (x,y) op(Marginal1(x),Marginal2(y));

Z3 = Joint(X,Y);
%figure;
%meshc(X,Y,Z3);


%% Generated correlated samples of arbitrary marginals using normal coupla
% family. Applies to higher dimensions


corr = 0.8;             %Select a Pearson correlation

% If spearmans correlation is give transformation is:
% corr = 2*sin(pi*spear/6);

%If Kendall is given:
% corr = sin(pi*kend/2);

a1 = rand(1,2000);       %Independant standard uniform
a2 = rand(1,2000);

w1 = norminv(a1);   %Transform to standard normal space
w2 = norminv(a2);   %Inverse cdf method

z1 = w1;
z2 = corr*w1+sqrt(1-corr^2)*w2; %Induce correlation in standard normal space using 

figure;
plot(z1,z2,'*');

figure;
plot(normcdf(z1),normcdf(z2),'*');

samples1 = logninv(normcdf(z1), 5, 1);
samples2 = betainv(normcdf(z2), 10,2);

%figure
%plot(z1,z2 ,'*');
figure
plot(samples1,samples2,'*');

%% Generating dependent samples using copula
% For this we must use conditioning sampling on the copula. First generate
% a random number, and find the conditional distribution of the copula at
% that point. This is the partial derrivative. Then sample this using it's
% inverse. The two sampled points will be distributed accoring to the
% copula. Transform to whichever marginal using it's inverse and wlala. 

Nsamples =2000;

a1 = rand(1,Nsamples);       %Independant standard uniform
a2 = rand(1,Nsamples);
a3 = zeros(1,Nsamples);


e = 0.0001;                 % Chose small number for numerical integration
spacing = 0.001;            % cdf descritization
numpoints=[0:spacing:1];
cdf = zeros(1,length(numpoints));

for i=1:Nsamples
    
    % Find marginal cdf
    %for l=1:length(numpoints)
    %   cdf(l)=(indep(numpoints(l),a1(i)+e)-indep(numpoints(l),a1(i)))./e;
    %end
    cdf = (C(numpoints,a1(i)+e)-C(numpoints,a1(i)))./e;
    %cdf = Gau(numpoints,a1(i),0.8);
    val = a2(i);
    iterand = -1;
    for k=1:length(cdf)     % Find location of sample in conditional
        if val<cdf(k)
            break
        end
    end
    a3(i)=numpoints(k);
end
figure
plot(a1,a3,'*')


% Transform samples to whichever marginal of choice


samples1 = norminv(a1);
samples2 = norminv(a3);

figure
plot(samples1,samples2,'*')

%figure
%plot(numpoints,cdf);

sigma2 = [1,-0.9;-0.9,1];

pp =mvnrnd(mu,sigma2,Nsamples);
figure
plot(pp(:,1),pp(:,2),'*');


figure;
ksdensity(pp(:,1))

figure;
ksdensity(pp(:,2))

figure
ksdensity([pp(:,1),pp(:,2)])



%% Choleski decomp method
a1 = randn(1,2000);
a2 = randn(1,2000);

figure
plot(a1,a2,'*');

A = [1,0.8;0.8,1];
L = chol(A,'lower');

Samples1 = L(1,1)*a1;
Samples2 = L(2,1)*a1 + L(2,2)*a2;

figure
plot(Samples1,Samples2,'*');

Samples1 = 600+Samples1;
Samples2 = 40+Samples2;

figure
plot(Samples1,Samples2,'*');




