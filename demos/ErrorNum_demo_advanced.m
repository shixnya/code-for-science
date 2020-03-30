%% Create two variables with covariance
x = 4 + randn(100, 1);
y = x + 0.5 * randn(100, 1);

%% Plot these two variables
figure;plot(x, y, '.');
axis([0, 8, 0, 8]);
setsize(2, 2);
box off
xlabel('x');
ylabel('y');
% print -dpng CovFig.png % print if necessary



%% calculate addition and subtraction without covariance
x_en = ErrorNum(mean(x), std(x));
y_en = ErrorNum(mean(x), std(y));
x_en + y_en
x_en - y_en


%% covariance calculation
cov_xy = cov(x, y);
x_en.covadd(y_en, cov_xy(1, 2))
x_en.covsub(y_en, cov_xy(1, 2))

%% ANOVA significance test
a = ErrorNum(3, 1);
b = ErrorNum(2, 1);
a.sig % this works
sig(a) % this also works
sig(a - b) % this also works
%(a - b).sig % this doesn't work.

%% Conversion
x = ErrorNum([1, 2, 3], 0.1 * ones(1, 3))

y = x.convert(@log, @(t) 1./t)

