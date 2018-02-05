clear; clc; close all;

img_dir = '~/Videos/director/';
images = dir([img_dir, '*.tiff']); 

fname = 'SP2_CL.avi';

v = VideoWriter(['~/Dropbox/', fname]);
open(v); 

for i = 1:numel(images)
    img = imread([img_dir, images(i).name]); 
    writeVideo(v, img)
end

close(v); 
