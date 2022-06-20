function [H,D,grid_h,grid_d,numVars]=compareData(dir_had,fname_dendro,timestep,stage,precision)

%fnameD=fname_dendro+"_"+string(timestep)+"_stage_"+string(stage)+".bin"
fnameD=fname_dendro+"_"+string(3*timestep+(stage+1))+".bin";
[H,grid_h,numVars]=xloadHAD(dir_had,timestep,stage,precision);
[D,grid_d,numVars]=xloaddata(fnameD,precision);

end