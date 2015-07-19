function [] = dialogDemo()
x = inputdlg('Enter space-separated numbers:',...
             'Enter Matrix', [5 50]);
data = str2num(x{:})';
disp(data);
disp(sum(data(:)));
disp(data');
end

