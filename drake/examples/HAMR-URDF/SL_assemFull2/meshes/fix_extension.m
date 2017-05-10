clear; clc; close all; 

files = dir('*.STL'); 

for i = 1:numel(files)
    [~, name, ext] = fileparts(files(i).name)
    oldname = [name, ext]; 
    newname = [name, '.stl'];
    movefile(oldname, newname); 
end
