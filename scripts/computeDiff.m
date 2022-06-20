[H_out,D,grid_h,grid_d,numVars]=compareData("../build/rhsvarsko","../build/rk3_step",0,0,'double');
E=H_out-D;
format long;
for i=1:24
    max(max(max(E(i,:,:,:))))
    %sprintf("var\t %d diff\t %f",i,max(max(max(E(i,:,:,:)))))
end