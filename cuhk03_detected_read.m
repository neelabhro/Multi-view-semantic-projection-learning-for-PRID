%CUHK03
%Neelabhro Roy
%IIIT-Delhi
                                                                                                                                                            
clear;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
clc;                                                                                                            
close all;
        
[srcdir,cam_pairs,num_pair] = folderList(fullfile('CUHK03/detected'));
camID = [];
personID = [];
imgs = {};        
cnt = 1;
        for cpair = 1:3 
            [~,cam_sets,num_cam] = folderList(fullfile(srcdir,cam_pairs{cpair}));
            for c = 1:num_cam
                [~,person_sets,num_person] = folderList(fullfile(srcdir,cam_pairs{cpair},cam_sets{c}));
                for p = 1:num_person
                    iminfo = dir(fullfile(srcdir,cam_pairs{cpair},cam_sets{c},person_sets{p},'*.png'));
                    for i = 1:numel(iminfo)
                        tmpim = imread(fullfile(srcdir,cam_pairs{cpair},cam_sets{c},person_sets{p},...
                            iminfo(i).name));
                        imgs{cnt} = imresize(tmpim,[128 64]);
                        camID(cnt) = str2num(cam_sets{c}(end))+10*cpair;
                        personID(cnt) = str2num(person_sets{p})+1000*cpair;
                        cnt = cnt + 1;
                    end
                end
            end
        end
