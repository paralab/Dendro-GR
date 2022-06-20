function [A,dimension,numVars] = xloaddata(filename,precision)
% [A, count] = xloaddata(filename, size, precision) loads
% a Matlab entity tfrom a binary file with certain precision.
% A can be a vector, a matrix, or a 3D volume
% size should be [M N] for M rows and N columns matrix
% precision controls the form and size of the result
% See the list of allowed precision under FREAD

fid = fopen(filename, 'rb');
dimension= [fread(fid,1,'int'),fread(fid,1,'int'),fread(fid,1,'int')];
numVars= fread(fid,1,'int');

A=zeros(numVars,dimension(1),dimension(2),dimension(3));

for v=1:numVars
     B=fread(fid,prod(dimension),precision);
     %size(B)
     for k=1:dimension(3)
         for j=1:dimension(2)
             for i=1:dimension(1)
                 A(v,i,j,k)=B((k-1)*dimension(2)*dimension(1)+(j-1)*dimension(1)+i,1);
             end    
         end
     end
end

fclose(fid);

end 