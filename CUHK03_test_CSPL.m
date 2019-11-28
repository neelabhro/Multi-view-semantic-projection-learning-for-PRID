%CUHK01 + LOMO
%Neelabhro Roy
%IIIT-Delhi
                                                                                                                                                            
clear;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
clc;                                                                                                            
close all;

feaFile = 'cuhk03-lomo-14670.mat';
%pcaFile = 'cuhk01100.mat';
pcaFile = 'cuhk03_pca_14670.mat';
numClass = 1367;
numFolds = 10;
numRanks = 100                                      ;                                                                     

%% load the extracted LOMO features
load(feaFile, 'descriptors');

p = randperm(1467);

camA1 = descriptors(:, 1:10:end);
%camA1 = camA1(:,1:1474);
camA2 = descriptors(:, 2:10:end);
camA3 = descriptors(:, 3:10:end);
camA4 = descriptors(:, 4:10:end);
camA5 = descriptors(:, 5:10:end);

camB1 = descriptors(:, 6:10:end);
camB2 = descriptors(:, 7:10:end);
camB3 = descriptors(:, 8:10:end);
camB4 = descriptors(:, 9:10:end);
camB5 = descriptors(:, 10:10:end);


A = 1:1:13670;

FirstElements  = A(:, 1:10:end);
SecondElements = A(:, 2:10:end);
ThirdElements  = A(:, 3:10:end);
FourthElements = A(:, 4:10:end);
FigthElements  = A(:, 5:10:end);

AllElementsMatrix = [FirstElements; SecondElements; ThirdElements; FourthElements; FigthElements];
[Rows, Cols] = size(AllElementsMatrix);
Array = reshape(AllElementsMatrix, [1, Rows*Cols]);
galFea1 = descriptors(:, Array);

FirstElements  = A(:, 6:10:end);
SecondElements = A(:, 7:10:end);
ThirdElements  = A(:, 8:10:end);
FourthElements = A(:, 9:10:end);
FigthElements  = A(:, 10:10:end);

AllElementsMatrix = [FirstElements; SecondElements; ThirdElements; FourthElements; FigthElements];
[Rows, Cols] = size(AllElementsMatrix);
Array = reshape(AllElementsMatrix, [1, Rows*Cols]);
probeFea1 = descriptors(:, Array);

%load(pcaFile, 'p');

%galFea1 = camA1(:, p(1:485));
gal1 = camA1(:, p(1368:1467));
gal2 = camA2(:, p(1368:1467));
gal3 = camA3(:, p(1368:1467));
gal4 = camA4(:, p(1368:1467));
gal5 = camA5(:, p(1368:1467));
galFea2 = (gal1 + gal2 + gal3 + gal4 + gal5)./5;

%probeFea1 = camB1(:, p(1:485));
probe1 = camB1(:, p(1368:1467));
probe2 = camB2(:, p(1368:1467));
probe3 = camB3(:, p(1368:1467));
probe4 = camB4(:, p(1368:1467));
probe5 = camB5(:, p(1368:1467));
probeFea2 = (probe1 + probe2 + probe3 + probe4 + probe5)./5;

%X1 = probeFea1;
%X2 = galFea1;
%X12 = probeFea2;
%X22 = galFea2;

    L = 0.000000000001;

    Lu = 0.05;    
    Lv = 0.2;
    La = 0.2;
    Lp = 0.2;
    Lw = 0.5;

    nu = 1;
    beta = 1;

    n = 6835;
    d = 100;
    k = d;
    
    TrainSet = zeros(13670,26960);
    TrainSet(1:6835 ,:) = galFea1';
    TrainSet(6836: end,:) = probeFea1';
    
    TestSet = zeros(200,26960);
    TestSet(1:100 ,:) = galFea2';
    TestSet(101: end,:) = probeFea2';
    
    %[X ,W] = matlabPCA(TrainSet',100);
    %TestSet(1:485 ,:) = galFea4';
    %TestSet(486: end,:) = probeFea4';
    
    load(pcaFile, 'X');
    load(pcaFile, 'W');
    
    X2 = X(:, 1:6835);
    X1 = X(:, 6836:end);
    %X1 = pca(probeFea1');
    %X2 = pca(galFea1');
    
    TestPCA = W' * TestSet';
    X22 = TestPCA(:, 1:100);
    X12 = TestPCA(:, 101:end);

    
    

    %X12 = pca(probeFea2');
    %X22 = pca(galFea2');

    % n is the size of the sample set
    % d is the feature dimension equal to 100
    % Components of Probe set X1 and gallery set X2 are matched image descriptors
    % K is the number of latent factors and is equal to d
    % U is the basis matrix capturing latent intrinsic structure of the input
    % data matrix
    % V1 and V2 (kxn) indicate the semantic representation of X1 and X2

    % X(dxn) = U(dxk)*V(kxn)
    %X1 = randi([0, 1], [d,n]);
    %X2 = randi([0, 1], [d,n]);

    U  = randi([0, 1], [d,k]);
    V1 = randi([0, 1], [k,n]);
    V2 = randi([0, 1], [k,n]);
    A  = randi([0, 1], [k,k]);
    P1 = randi([0, 1], [k,d]);
    P2 = randi([0, 1], [k,d]);
    W2 = eye(k);
    



%% Main algorithm
    for i = 1:500
        U  = (( W2*X1 * transpose(V1)) + ( W2*X2 * transpose(V2)))/((( V1 * transpose(V1)) + ( V2 * transpose(V2)) + (Lu*eye(k))));
        V1 = (((transpose(U) * U) + (nu + beta + Lv) * eye(k))) \ ((transpose(U) *W2* X1) + (beta* A * V2) + nu * P1 * W2*X1);
        V2 = (((transpose(U) * U) + ( beta * transpose(A) * A) + (nu + Lv) .* eye(k))) \ ((transpose(U) *W2* X2) + (beta* transpose(A) * V1) + nu * P2*W2 * X2);
        P1 = (V1 * transpose(X1)) / ((X1 * transpose(X1)) + (Lp/nu)*eye(k));
        P2 = (V2 * transpose(X2)) /((X2 * transpose(X2)) + (Lp/nu)*eye(k));
        A  = (V1 * transpose(V2)) /((V2 * transpose(V2)) + (La/beta)*eye(k));
        %W2 = inv((2*X1*X1' + 2*X2*X2' + Lw*eye(k)))*(X1*V1'*U' + X2*V2'*U' + X1*V1'*P1 + X2*V2'*P2);
        %W2 = inv((4*X1*X1' + 4*X2*X2' + Lw*eye(k) -X1*X2' -X2*X1'))*(X1*V1'*U' + X2*V2'*U' + X1*V1'*P1 + X2*V2'*P2);
        %W2 = W2';
        %W2 = eye(k);
        %for j = 1:50
            
        %    Vw = 4*W2'*(X1*X1') + 4*W2'*(X2*X2') - 2*(X1*V1')*U'  - 2*(X2*V2')*U' - 2*nu*X1*V1'*P1' - 2*nu*X2*V2'*P2' + 2*nu*(P1'*P1)*W2'*(X1*X1') + 2*nu*(P2'*P2)*W2'*(X2*X2') + 2*Lw*W2' - X1*X2' - X2*X1';
        %    W2 = W2 - L*Vw;
            
            %W=W./norm(W);
        %end 
        %W2 = W2';
        
    end

    
    D = zeros(100,100);
    for m = 1:100
    
        v1 = P1*W2*(X12(:,m));
    
        for i = 1:100
            v2 = P2*W2*(X22(:,i));
            D(m,i) = norm(((v1 - A*v2)));
        end
        
    end
    
CMC(D,100);