function display_diff(H, D, var, slice,dir)
% function display_diff(H, D, var, slice)

E = H - D;
if(dir==3)
    subplot(2,2,1), imagesc(squeeze(E(var, :,:,slice))); colorbar; title('Error'); % [-1e-3, 1e-3]
    subplot(2,2,2), imagesc(squeeze(H(var, :,:,slice))); colorbar; title('HAD');
    subplot(2,2,3), imagesc(squeeze(D(var,:,:,slice))); colorbar; title('Dendro'); 
elseif(dir==2)
    subplot(2,2,1), imagesc(squeeze(E(var, :,slice,:))); colorbar; title('Error'); % [-1e-3, 1e-3]
    subplot(2,2,2), imagesc(squeeze(H(var, :,slice,:))); colorbar; title('HAD');
    subplot(2,2,3), imagesc(squeeze(D(var,:,slice,:))); colorbar; title('Dendro'); 
elseif(dir==1)
    
    subplot(2,2,1), imagesc(squeeze(E(var, slice,:,:))); colorbar; title('Error'); % [-1e-3, 1e-3]
    subplot(2,2,2), imagesc(squeeze(H(var, slice,:,:))); colorbar; title('HAD');
    subplot(2,2,3), imagesc(squeeze(D(var,slice,:,:))); colorbar; title('Dendro'); 
    
end


clear E;
end