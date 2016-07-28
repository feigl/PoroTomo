

clear all; close all; clc

File_url{1} = 'https://gdr.openei.org/files/828/Brady_Obs_Metadata.csv';
File_url{2} = 'https://gdr.openei.org/files/827/Reftek_metadata.csv';
File_url{3} = 'https://gdr.openei.org/files/826/Nodal_continuous_metadata%20(1).csv';

File_names{1} = {'Brady_Obs_Metadata.csv'};
File_names{2} = {'Reftek_metadata.csv'};
File_names{3} = {'Nodal_continuous_metadata%20(1).csv'};

for i = 1:numel(File_url)
    urlwrite(File_url{i},char(File_names{i}));
end


for i = 1:numel(File_url)
    csv2struct(char(File_names{i}))
end


    






