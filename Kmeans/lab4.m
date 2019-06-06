clc
clear all
close all
%% Set Up

I= imread('house.tiff');
X = reshape(I, 256*256, 3);
X = double(X);
R = [0,0,1];G = [0,1,0];B = [1,0,0];
figure, plot3(X(:,1), X(:,2), X(:,3),'.','Color',[1, 0, 0])

mean1 = [.1 .1 .1];

%% a- c=2

inimean = [.5 .6 .7
     .8 .8 .9];

inimean = inimean*256;
updatemean = zeros(size(inimean));
J = [];
m1 = [inimean(1,:)];
m2 = [inimean(2,:)];

while (updatemean ~= inimean)

    updatemean = inimean;
    
    j1 = (X - repmat(inimean(1,:), size(X,1), 1));
    j2 = (X - repmat(inimean(2,:), size(X,1), 1));
    j1 = sum(j1.^2, 2);    
    j2 = sum(j2.^2, 2);

    cluster1 = j1 < j2;
    cluster2 = ~cluster1;

    inimean(1,:) = sum(X(cluster1, :)) / sum(cluster1);
    inimean(2,:) = sum(X(cluster2, :)) / sum(cluster2);
    J = [J sum(min(j1, j2))];
    m1 = [m1; inimean(1,:)];
    m2 = [m2; inimean(2,:)];
    
  
    
    
    
end
inimean
m1
m2
%%

close all
plot(J)
title('Error Criterion');

figure
plot3(m1(:, 1), m1(:, 2), m1(:, 3), '-o')
hold all
plot3(m2(:, 1), m2(:, 2), m2(:, 3), '-o')
title('Clustering Means Two Stages')
%%
figure
x0 = X(cluster1, :);
x00 = X(cluster2, :);
plot3(x0(:,1), x0(:,2), x0(:,3),'+','Color', inimean(1,:)/256)
hold all
plot3(x00(:,1), x00(:,2), x00(:,3),'+','Color', inimean(2,:)/256)
title('Data Samples RGB Space');
%%

figure
x1 = repmat(inimean(1,:), size(X,1), 1) .* repmat(cluster1, 1, size(X,2));
x1 = x1 + repmat(inimean(2,:), size(X,1), 1) .* repmat(cluster2, 1, size(X,2));
x1 = reshape(x1, size(I, 1), size(I, 2), 3);
subplot(1,2,1);imshow(I)
title('Original');
subplot(1,2,2);imshow(x1/256)
title('Labeled Form');

%% b

c = 5;

inib = [
  161.25    61.50    171.25
  151.25    30.50    161.25
  240.50    170.75   100.25
  45.50     120.0    200.65
  170.25    110.75   95.75
];
XX=[];
meanbb = inib;

upm = zeros(size(meanbb));

while (upm ~= meanbb)

    upm = meanbb;
    
    J = zeros(size(X,1), c);
    
    for i = [1:c]
        j1 = (X - repmat(meanbb(i,:), size(X,1), 1));
        j1 = sum(j1.^2, 2);
        J(:,i) = j1;
    end
    [a, cluster] = min(J, [], 2);
    
    for i = [1:c]
        tc = (cluster==i);
        meanbb(i, :) = sum(X(tc, :)) / sum(tc);

    end
XX = zeros(size(X));
for i = [1:c]
    tc = (cluster==i);
    XX = XX + repmat(meanbb(i,:), size(X,1), 1) .* repmat(tc, 1, size(X,2));
end
XX = reshape(XX, size(I, 1), size(I, 2), 3);
subplot(1,2,1);imshow(I)
subplot(1,2,2);imshow(XX/256)



end

mbb = meanbb
cluster1 = cluster;
figure
for i = [1:c]
    tc = (cluster==i);
    Xi = X(tc, :);
    plot3(Xi(:,1), Xi(:,2), Xi(:,3),'.','Color', meanbb(i,:)/256)
    hold all
end
title('RGB Cluster');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
inib2 = [
  200.50  250.75  165.25
  190.5   225.70  45.25
  180.5   245.25  250.35
  60.5   135.25   110.75
  150.25  30.5    220.75
];

meanbb2 = inib2;

upm = zeros(size(meanbb2));

while (upm ~= meanbb2)

    upm = meanbb2;
    
    J = zeros(size(X,1), c);
    
    for i = [1:c]
        j2 = (X - repmat(meanbb2(i,:), size(X,1), 1));
        j2 = sum(j2.^2, 2);
        J(:,i) = j2;
    end
    [a, cluster] = min(J, [], 2);
    
    for i = [1:c]
        tc = (cluster==i);
        meanbb2(i, :) = sum(X(tc, :)) / sum(tc);
    end
    
end
mb2 = meanbb2
cluster2 = cluster;

%
figure
for i = [1:c]
    tc = (cluster==i);
    Xi = X(tc, :);
    plot3(Xi(:,1), Xi(:,2), Xi(:,3),'.','Color', meanbb2(i,:)/256)
    hold all
end
title('RGB Cluster 2');
%% c

N = size(X,1);
q1 = 0;

for i = [1:c]
    tc = (cluster1==i);
    Xi = X(tc, :);
    mj = sort(sum((mbb - repmat(mbb(i,:), c, 1)).^2, 2).^.5);
    q1 = q1 + sum(sum((Xi - repmat(mbb(i,:), size(Xi,1), 1)).^2, 2).^.5) / mj(2);
end

q1 = q1 / N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
q2 = 0;

for i = [1:c]
    tc = (cluster2==i);
    Xi = X(tc, :);
    mj = sort(sum((mb2 - repmat(mb2(i,:), c, 1)).^2, 2).^.5);
    q2 = q2 + sum(sum((Xi - repmat(mb2(i,:), size(Xi,1), 1)).^2, 2).^.5) / mj(2);
end

q2 = q2 / N