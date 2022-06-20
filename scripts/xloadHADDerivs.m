function [A,dimension,numVars] =xloadHADDerivs(filename,timestep,stage,precision)
%dName dir name
%numsteps 
% precision

% enum VAR {U_ALPHA=0,U_CHI,U_K,U_GT0,U_GT1,U_GT2,U_BETA0,U_BETA1,U_BETA2,U_B0,U_B1,U_B2,U_SYMGT0,U_SYMGT1,U_SYMGT2,U_SYMGT3,U_SYMGT4,U_SYMGT5,U_SYMAT0,U_SYMAT1,U_SYMAT2,U_SYMAT3,U_SYMAT4,U_SYMAT5};
%fnames=["dt_alpha","dt_chi","dt_trK","dt_Gamt1","dt_Gamt2","dt_Gamt3","dt_shift1", "dt_shift2","dt_shift3", "dt_gb1","dt_gb2","dt_gb3", "dt_gt11","dt_gt12","dt_gt13", "dt_gt22", "dt_gt23","dt_gt33","dt_A11","dt_A12", "dt_A13", "dt_A22","dt_A23","dt_A33"];
%fnames=["alpha","chi","trK","Gamt1","Gamt2","Gamt3","shift1", "shift2","shift3", "gb1","gb2","gb3", "gt11","gt12","gt13", "gt22", "gt23","gt33","A11","A12", "A13", "A22","A23","A33"];
%file_index=3*timestep+stage;
numVars=1;
fid = fopen(filename, 'rb');
%q = fread(fid,1,'uint')
x = fread(fid,1,'uint');
y = fread(fid,1,'uint');
z = fread(fid,1,'uint');
dimension= [x, y, z];
fclose(fid);


A=zeros(numVars,dimension(1),dimension(2),dimension(3));
    
for v=1:numVars
    fid = fopen(filename, 'rb');
    x = fread(fid,1,'uint');
    y = fread(fid,1,'uint');
    z = fread(fid,1,'uint');
    B=fread(fid,prod(dimension),precision);
     for k=1:dimension(3)
         for j=1:dimension(2)
             for i=1:dimension(1)
                 A(v,i,j,k)=B((k-1)*dimension(2)*dimension(1)+(j-1)*dimension(1)+i,1);
             end    
         end
     end
     fclose(fid);
end
    
end
