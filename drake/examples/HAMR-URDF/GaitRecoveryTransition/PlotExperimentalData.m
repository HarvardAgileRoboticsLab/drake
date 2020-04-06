clear; clc; close all;

datadir = '~/Dropbox/GaitRecoveryandTransition/data/ExpData/HypII';
files = dir([datadir '/**/*.mat']); %dir(fullfile(datadir, '*_30Hz_*.mat'));

data = cell(numel(files), 1);
for i = 1:numel(files)
    fname = [files(i).folder, '/', files(i).name];
    data{i} = load(fname);
    tok = strsplit(files(i).name, '_');
    freq(i) = str2double(tok{7}(1:end-2));
end

%%
DC = cellfun(@(x) x.params.dutycycle, data);
DL = cellfun(@(x) x.params.push, data);
% freq = cellfun(@(x) x.params.freq, data);


[freq_sort, ind_sort] = sort(freq, 'ascend');
DC_sort = DC(ind_sort);
DL_sort = DL(ind_sort);
data_sort = data(ind_sort);

Nf = numel(unique(freq_sort));
Ndc = numel(unique(DC_sort));
Ndl = numel(unique(DL_sort));

freq_mat = reshape(freq_sort, Ndc, Ndl, Nf);
DC_mat = permute(reshape(DC_sort, Ndl, Ndc, Nf), [2, 1, 3]);
DL_mat = permute(reshape(DL_sort, Ndl, Ndc, Nf), [2, 1, 3]);
data_mat = permute(reshape(data_sort, Ndl, Ndc, Nf), [2, 1, 3]);

avg_speed = zeros(Ndc, Ndl, Nf);
rms_err =  zeros(Ndc, Ndl, Nf);

[b, a] = butter(2, 0.05, 'low'); 
for i = 1:Nf % frequency
    for j = 1:Ndc % swing duty cycle (rows)
        for k = 1:Ndl % lift push (columns)
            
            % unpack parameters
            datai = data_mat{j, k, i};
            RAMPCYC = datai.params.initcyc;
            DATACYC = datai.params.datacyc;
            
%             fprintf('(%d, %d) is %d swing and %d lift \r', j, k, ...
% %                 DC_mat(j,k),  DL_mat(j,k))
%             disp(datai.params.dutycycle)
%             disp(datai.params.push)
%             disp(datai.params.freq)
            
            tt_sol = datai.savdata.t;
            dt_sol = mean(diff(tt_sol));
            xyz_sol = datai.savdata.xyz;
                 
            figure(i); clf; hold on; 
            plot(datai.savdata.xleg_vic(:,[1,3])); hold on; 
            title(sprintf('DL: %f, DC %f', datai.params.push, datai.params.dutycycle))
            % compute avg speed and error
            ind0 = find(tt_sol > RAMPCYC/freq_mat(j,k,i), 1, 'first');
            ind1 = find(tt_sol >= (RAMPCYC+DATACYC)/freq_mat(j,k,i), 1, 'first');
            if isempty(ind1); ind1 = numel(tt_sol); end
            x = xyz_sol(ind0:ind1,1) - xyz_sol(ind0,1);
%             ind_dropped = find(x < 0)
%             x(ind_dropped) = x(ind_dropped-1);
            xf = filtfilt(b, a, x); 
            y = xyz_sol(ind0:ind1,2) - xyz_sol(ind0,2);
%             avg_speed(j,k,i) = sum(sqrt(diff(x).^2 + diff(y).^2))/(tt_sol(end) - tt_sol(ind0));
            avg_speed(j,k,i) = sum(x(end)-x(1))/(tt_sol(end) - tt_sol(ind0));
            rms_err(j,k,i) = rms(datai.savdata.pos_err(:));
            
%             figure(i); hold on;
%             title([num2str(freq_mat(1,1,i)) ' Hz'])
%             plot(tt_sol(ind0:ind1), x(:, 1));
        end
    end
end

for i = 1:Nf
    figure(Nf+i); clf;
    subplot(2,1,1); hold on;
    title([num2str(freq_mat(1,1,i)) ' Hz'])
    contourf(DC_mat(:,:,i), DL_mat(:,:,i), avg_speed(:,:,i)./freq_mat(:,:,i))
    xlabel('swing dc')
    ylabel('lift dc')
    colorbar
    
    subplot(2,1,2); hold on;
    title([num2str(freq_mat(1,1,i)) ' Hz'])
    contourf(DC_mat(:,:,i), DL_mat(:,:,i), rms_err(:,:,i))
    colorbar    
    
end

tilefigs; 




