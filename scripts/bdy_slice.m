close all;
figure

for i=1:65
    imagesc(squeeze(D(1,i,:,:))), colorbar,title(iD);
    pause(1);
end 