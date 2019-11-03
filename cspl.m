%Neelabhro Roy
%IIIT, Delhi

clear;
clc;
close all;

imgDir = 'CUHK03/detected/';
addpath('../bin/');
%% 1_001_1_01.png means the 1st set, the 001th person, the 1st camera, and the 01th image.

%% Image list
%list = dir('CUHK03/detectedlist = dir('CUHK03/detected/1*.png'));
list = dir(['CUHK03/detected/*.png']);
n = length(list);

%% Read images
%for i = 1 : n
%    info = imfinfo([imgDir, list(i).name]);
%    images = zeros(info.Height, info.Width, 3, n, 'uint8');
%    images(:,:,:,i) = imread([imgDir, list(i).name]);
%    images1(:,:,:,i) = imresize(images(:,:,:,i),[128 48]);
%end

%imshow(images1(:,:,:,1));

    Lu = 0.05*1;
    L = 0.000000001;
    Lv = 0.2*1;
    La = 0.2*1;
    Lp = 0.2*1;
    Lw = 05*1;
    nu = 1;
    beta = 1;
    
    
    n = 100;
    %843 Persons from Set 1
    %842 Persons for simplicity
    %1 image from each of the 2 views to make probe and gallery
    d = 100;
    k = d;

    X1 = randi([0, 1], [d,n]);
    X2 = randi([0, 1], [d,n]);
    X3 = randi([0, 1], [d,n]);
    X4 = randi([0, 1], [d,n]);
    X5 = randi([0, 1], [d,n]);
    X6 = randi([0, 1], [d,n]);
    
    
    U  = randi([0, 1], [d,k]);
    V1 = randi([0, 1], [k,n]);
    V2 = randi([0, 1], [k,n]);
    V3 = randi([0, 1], [k,n]);
    V4 = randi([0, 1], [k,n]);
    V5 = randi([0, 1], [k,n]);
    V6 = randi([0, 1], [k,n]);    
    
    A12  = randi([0, 1], [k,k]);
    A13  = randi([0, 1], [k,k]);
    A14  = randi([0, 1], [k,k]);
    A15  = randi([0, 1], [k,k]);
    A16  = randi([0, 1], [k,k]);

    A23  = randi([0, 1], [k,k]);
    A24  = randi([0, 1], [k,k]);
    A25  = randi([0, 1], [k,k]);
    A26  = randi([0, 1], [k,k]);

    A34  = randi([0, 1], [k,k]);    
    A35  = randi([0, 1], [k,k]);
    A36  = randi([0, 1], [k,k]);
    
    A45  = randi([0, 1], [k,k]);
    A46  = randi([0, 1], [k,k]);
    
    A56  = randi([0, 1], [k,k]);
    
    P1 = randi([0, 1], [k,d]);
    P2 = randi([0, 1], [k,d]);
    P3 = randi([0, 1], [k,d]);
    P4 = randi([0, 1], [k,d]);
    P5 = randi([0, 1], [k,d]);
    P6 = randi([0, 1], [k,d]);
    W = eye(k);


    
%% Main algorithm
    for i = 1:50
        U  = (( W'*X1 * transpose(V1)) + ( W'*X2 * transpose(V2)))/((( V1 * transpose(V1)) + ( V2 * transpose(V2)) + (Lu*eye(k))));

        V1 = (((transpose(U) * U) + (nu + 5*beta + Lv) * eye(k))) \ ((transpose(U) *W'* X1) + (beta* A12 * V2) + (beta* A13 * V3) + (beta* A14 * V4) + (beta* A15 * V5) + (beta* A16 * V6) + nu * P1 * W'*X1);
        V2 = (((transpose(U) * U) + ( beta * transpose(A12) * A12) + (nu + Lv + 4*beta) .* eye(k))) \ ((transpose(U) *W'* X2) + (beta* transpose(A12) * V1) + (beta* A23* V3) + (beta* A24* V4) + (beta* A25* V5) + (beta* A26* V6) + nu * P2*W'* X2);
        V3 = (((transpose(U) * U) + ( beta * transpose(A23) * A23) + ( beta * transpose(A13) * A13) + (nu + Lv + 3*beta) .* eye(k))) \ ((transpose(U) *W'* X3) + (beta* transpose(A23) * V2) + (beta* A34* V4) + (beta* A35* V5) + (beta* A36* V6) + nu * P3'*W'* X3);
        V4 = (((transpose(U) * U) + ( beta * transpose(A34) * A34) + ( beta * transpose(A14) * A14) + ( beta * transpose(A24) * A24) + (nu + Lv + 2*beta) .* eye(k))) \ ((transpose(U) *W'* X4) + (beta* transpose(A34) * V3) + (beta* A45* V5) + (beta* A46* V6) + nu * P4'*W'* X4);
        V5 = (((transpose(U) * U) + ( beta * transpose(A45) * A45) + ( beta * transpose(A15) * A15) + ( beta * transpose(A25) * A25) + ( beta * transpose(A35) * A35) + (nu + Lv + 1*beta) .* eye(k))) \ ((transpose(U) *W'* X5) + (beta* transpose(A45) * V4) + (beta* A56* V6) + nu * P5'*W'* X5);
        V6 = (((transpose(U) * U) + ( beta * transpose(A16) * A16) + ( beta * transpose(A26) * A26) + ( beta * transpose(A36) * A36) + ( beta * transpose(A46) * A46) + ( beta * transpose(A56) * A56) + (nu + Lv + 0*beta) .* eye(k))) \ ((transpose(U) *W'* X6) + (beta* transpose(A16) * V1) + (beta* A26'* V2) + (beta* A36'* V3) + (beta* A46'* V4) + (beta* A56'* V5) + nu * P6'*W'* X6);
        
        P1 = (V1 * transpose(X1)) / ((X1 * transpose(X1)) + (Lp/nu)*eye(k));
        P2 = (V2 * transpose(X2)) / ((X2 * transpose(X2)) + (Lp/nu)*eye(k));
        P3 = (V3 * transpose(X3)) / ((X3 * transpose(X3)) + (Lp/nu)*eye(k));
        P4 = (V4 * transpose(X4)) / ((X4 * transpose(X4)) + (Lp/nu)*eye(k));
        P5 = (V5 * transpose(X5)) / ((X5 * transpose(X5)) + (Lp/nu)*eye(k));
        P6 = (V6 * transpose(X6)) / ((X6 * transpose(X6)) + (Lp/nu)*eye(k));
        
        A12  = -2*beta*V1*V2' + 2*A12*(V2)*V2' + 2*La*A12;
        A13  = -2*beta*V1*V3' + 2*A13*(V3)*V3' + 2*La*A13;
        A14  = -2*beta*V1*V4' + 2*A14*(V4)*V4' + 2*La*A14;
        A15  = -2*beta*V1*V5' + 2*A15*(V5)*V5' + 2*La*A15;
        A16  = -2*beta*V1*V6' + 2*A16*(V6)*V6' + 2*La*A16;
        A23  = -2*beta*V2*V3' + 2*A23*(V3)*V3' + 2*La*A23;
        A24  = -2*beta*V2*V4' + 2*A24*(V4)*V4' + 2*La*A24;
        A25  = -2*beta*V2*V5' + 2*A25*(V5)*V5' + 2*La*A25;
        A26  = -2*beta*V2*V6' + 2*A26*(V6)*V6' + 2*La*A26;
        A34  = -2*beta*V3*V4' + 2*A34*(V4)*V4' + 2*La*A34;
        A35  = -2*beta*V3*V5' + 2*A35*(V5)*V5' + 2*La*A35;
        A36  = -2*beta*V3*V6' + 2*A36*(V6)*V6' + 2*La*A36;
        A45  = -2*beta*V4*V5' + 2*A45*(V5)*V5' + 2*La*A45;
        A46  = -2*beta*V4*V6' + 2*A46*(V6)*V6' + 2*La*A46;
        A56  = -2*beta*V5*V6' + 2*A56*(V6)*V6' + 2*La*A56;
        
        for j = 1:5
            Vw = 2*(X1)*X1'*W + 2*(X2)*X2'*W + 2*(X3)*X3'*W + 2*(X4)*X4'*W + 2*(X5)*X5'*W + 2*(X6)*X6'*W - 2*(X1*V1')*U' - 2*(X2*V2')*U' - 2*(X3*V3')*U' - 2*(X4*V4')*U' - 2*(X5*V5')*U' - 2*(X6*V6')*U' - 2*nu*X1*V1'*P1' - 2*nu*X2*V2'*P2' - 2*nu*X3*V3'*P3' - 2*nu*X4*V4'*P4' - 2*nu*X5*V5'*P5' - 2*nu*X6*V6'*P6' + 2*nu*(X1)*X1'*W*(P1'*P1) + 2*nu*(X2)*X2'*W*(P2'*P2) + 2*nu*(X3)*X3'*W*(P3'*P3) + 2*nu*(X4)*X4'*W*(P4'*P4) + 2*nu*(X5)*X5'*W*(P5'*P5) + 2*nu*(X6)*X6'*W*(P6'*P6);
            W = W - L*Vw;
        end
    end

    
    D = zeros(n,n);
    for m = 1:n
    
        v1 = P1*W*(X1(:,m));
    
        for i = 1:n
            v2 = P2*W*(X2(:,i));
            D(m,i) = norm(((v1 - A12*v2)));
        end
        
    end
    
CMC(D,100);
