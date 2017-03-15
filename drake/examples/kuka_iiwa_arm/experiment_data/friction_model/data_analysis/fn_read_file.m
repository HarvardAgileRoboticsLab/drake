function ret = fn_read_file(path, file_name, num_col)
    fname = sprintf('%s/%s.txt', path, file_name);
    fileID = fopen(fname,'r');
    
    ret = zeros(num_col, 1);
    j = 1;
    
    while(~feof(fileID))
        for i = 1:num_col
            C = textscan(fileID, '%f',1);
            if(feof(fileID))
               break; 
            end
            ret(i,j) = C{1};
        end
        j = j+1;
    end
    fclose(fileID);
end