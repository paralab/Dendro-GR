close all;
figure

var=1;
slice=32;
T=0;
S=4;

for t=0:T
    
    for s=1:S
    
        [H,D,grid_h,grid_d,numVars]=compareData("../build/haddata","../build/rkU",t,s,'double');
        for var=1:24
            
            I = squeeze(H(var,:,:,slice));
            J = squeeze(D(var,:,:,slice));
        
            subplot(1, 3, 1), imagesc(I); title("HAD");
            subplot(1, 3, 2), imagesc(J); title("Dendro");colorbar;
            str=sprintf("var %d t: %d s: %d id: %d ",var,t,s,(3*t+(s+1)));
            subplot(1, 3, 3), imagesc((I-J)); title(str); colorbar;
            %set(gca,'colorscale','log')
        
            pause;
            
        end
    
        
    end
    
end

