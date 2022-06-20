dxAlpha=xloadHADDerivs("../build/derivs/adx_alpha_0.dat",0,0,'double');
dyAlpha=xloadHADDerivs("../build/derivs/ady_alpha_0.dat",0,0,'double');
dzAlpha=xloadHADDerivs("../build/derivs/adz_alpha_0.dat",0,0,'double');

dxAlpha_had=squeeze(dxAlpha(1,:,:,:));
dyAlpha_had=squeeze(dyAlpha(1,:,:,:));
dzAlpha_had=squeeze(dzAlpha(1,:,:,:));

[H_out,D,grid_h,grid_d,numVars]=compareData("../build/rhsvars","../build/rk3_step",0,0,'double');

dxAlpha_dendro=squeeze(D(7,:,:,:));
dyAlpha_dendro=squeeze(D(8,:,:,:));
dzAlpha_dendro=squeeze(D(9,:,:,:));

diff_x=(dxAlpha_dendro-dxAlpha_had);
diff_y=(dyAlpha_dendro-dyAlpha_had);
diff_z=(dzAlpha_dendro-dzAlpha_had);

% dx=20/64;
% x=-10;
% 
% %D_matlab=zeros([65,65,65]);
% i=1;
% for k=1:65
%     z=-10+(k-1)*dx;
%     for j=1:65
%         y=-10+(j-1)*dx;
%         inv_r=1.0/sqrt((x*x)+(y*y)+(z*z));
%         D(1,i,j,k)=-inv_r*(x*dxAlpha_had(i,j,k) + y* dyAlpha_had(i,j,k) + z * dzAlpha_had(i,j,k) + 1.0 *((1.0 - 0.25*sin((31.0/17.0)*x))-1.0) );
%     end
% end
% 
% 
% I=squeeze(D(1,1,:,:));
% J=squeeze(H_out(1,1,:,:));
% 
% [H_out,D,grid_h,grid_d,numVars]=compareData("../build/data","../build/rk3_step",0,0,'double');
% M=squeeze(D(1,1,:,:));
