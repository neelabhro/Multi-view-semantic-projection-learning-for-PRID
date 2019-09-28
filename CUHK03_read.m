clc; clear;

load('/home/kevinh/dataset/reid/cuhk03_release/cuhk-03.mat');

detected_or_labeled = 0;

if detected_or_labeled == 1
    img_vector = detected;
    load('cuhk03_multishot_config_detected.mat');
    path  = '/home/kevinh/dataset/reid/pytorch/CUHK-03/detected/';
else
    img_vector = labeled;
    load('cuhk03_multishot_config_labeled.mat');
    path  = '/home/kevinh/dataset/reid/pytorch/CUHK-03/labeled/';
end

for i = 1: length(filelist)
    img_name = filelist{i};
    img_name = strrep(img_name,'png','jpg');
    j1 = str2num(img_name(1));      % camera pair index
    j2 = str2num(img_name(3:5));    % identity index
    j3 = str2num(img_name(7));      % camera A or B of this pair
    j4 = str2num(img_name(9:10));   % image index of this identity
    if any(train_idx==i)
        train_all_path = fullfile(path,'train_all',num2str(labels(i)));
        if ~exist(train_all_path,'dir')
            mkdir(train_all_path);
        end
        imwrite(img_vector{j1}{j2,j4}, fullfile(train_all_path,img_name));
        save_path = fullfile(path,'train',num2str(labels(i)));
        if ~exist(save_path,'dir')
            val_path  = fullfile(path,'val',num2str(labels(i)));
            mkdir(val_path);
            imwrite(img_vector{j1}{j2,j4}, fullfile(val_path,img_name));
            mkdir(save_path);
            continue
        end
    elseif any(gallery_idx==i)
        save_path = fullfile(path,'gallery',num2str(labels(i)));
    elseif any(query_idx==i)
        save_path = fullfile(path,'query',num2str(labels(i)));
    end
    if ~exist(save_path,'dir')
        mkdir(save_path);
    end
    imwrite(img_vector{j1}{j2,j4}, fullfile(save_path,img_name));
    disp(i)
end
