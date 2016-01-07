A=rand(5,5);
B=rand(2,9);

fileID=fopen('test.bin','w');
fwrite(fileID, A, 'double');
fwrite(fileID, B, 'double');
fclose(fileID);

fileID=fopen('test.bin','r');
C=fread(fileID, size(A), 'double');
D=fread(fileID, size(B), 'double');
fclose(fileID);

disp(A-C)
disp(B-D)