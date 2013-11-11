close all
clear all 
clc

filename = 'overview.txt';

%[A,B,C] = reader(file)

%file = fopen(['data/',filename]);

%text = cellstr(file)

%Data = fscanf(file, '%c2',inf)

%Data = strjoin(Data, ',')

%Data = textscan(file, '%s,%s,%s,%s,%s,%s,%s,%s,%s')

file = fileread(['data/',filename]);

i = 1;
str = [];

while i < length(file)
   if file(i) == ','
       str = [str, ' '];
   elseif file(i) == ' '
       
   elseif file(i) == char(10)
       %disp('newline')
       str = [str, file(i)];
   else 
       str = [str, file(i)];
   end
 i = i+1   ;
end

[A,B,C,D,E,F,G,H,I] = strread(str,'%s %s %s %s %s %s %s %s %s');

N = length(A);

Data = NaN(N-1,8);
i =2
while i < N
    if strcmp(B{i},R) | strcmp(B{i},Rx)| strcmp(B{i},Ry)| strcmp(B{i},Rxy)
       Data(i-1,3) = 'Rsd';
    elseif strcmp(B{i},Angle) | strcmp(B{i},Anglexy) | strcmp(B{i},Angleyz) | strcmp(B{i},Anglexz)
       Data(i-1,3) = 'Asd';
    elseif strcmp(B{i},L) | strcmp(B{i},Lxy)| strcmp(B{i},Lyz) | strcmp(B{i},Lxz)
       Data(i-1,3) = 'Lsd';
    else
       Data(i-1,3) = 'Psd';
    end
i=i+1
end




for i = 2:N 
    Data(i-1,1) = strread([C{i}],'%f %*s');
    Data(i-1,2) = strread([G{i}],'%f %*s');
    
    Data(i-1,5) = strread([D{i}],'%f %*s');
    Data(i-1,6) = strread([E{i}],'%f %*s');
    Data(i-1,7) = strread([F{i}],'%f %*s');

end

