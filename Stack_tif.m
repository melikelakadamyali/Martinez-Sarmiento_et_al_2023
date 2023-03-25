
[data_load,path] = load_stack_tiff_file();

for i = 1:length(data_load)
    imwrite(data_load{i}.image, [path,'stacked_image.tif'],'writemode','append','Compression','none');
end


function [data_load,path] = load_stack_tiff_file()
[file_name,path] = uigetfile({'*.tif';'*.tiff'},'Select TIFF File(s)','MultiSelect','on');
if isequal(file_name,0)
    data_load=[];
else
    file_name=cellstr(file_name);
    for i=1:size(file_name,2)
        image = imread(fullfile(path,file_name{1,i}),1);
        data_load{i}.image = image;
        data_load{i}.name = file_name{1,i};
        data_load{i}.type = 'image';
        data_load{i}.info = 'NaN';        
        clear image
    end
end
end