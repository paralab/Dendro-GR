if(V_ID == 0) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_alpha + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_alpha + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_alpha + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_alpha + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_alpha + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_alpha + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_xx<pw,pencils,pencil_sz>(deriv_evars->grad2_0_0_alpha + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_yy<pw,pencils,pencil_sz>(deriv_evars->grad2_1_1_alpha + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_zz<pw,pencils,pencil_sz>(deriv_evars->grad2_2_2_alpha + Z_ID * BLK_SZ, su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_0_alpha + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad2_0_1_alpha + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_0_2_alpha + Z_ID * BLK_SZ , su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_1_alpha + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_1_2_alpha + Z_ID * BLK_SZ , su, blk);
return;
} // end of var alpha 
if(V_ID == 1) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_chi + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_chi + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_chi + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_chi + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_chi + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_chi + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_xx<pw,pencils,pencil_sz>(deriv_evars->grad2_0_0_chi + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_yy<pw,pencils,pencil_sz>(deriv_evars->grad2_1_1_chi + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_zz<pw,pencils,pencil_sz>(deriv_evars->grad2_2_2_chi + Z_ID * BLK_SZ, su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_0_chi + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad2_0_1_chi + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_0_2_chi + Z_ID * BLK_SZ , su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_1_chi + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_1_2_chi + Z_ID * BLK_SZ , su, blk);
return;
} // end of var chi 
if(V_ID == 2) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_K + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_K + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_K + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_K + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_K + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_K + Z_ID * BLK_SZ , su, blk);
return;
} // end of var K 
if(V_ID == 3) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_Gt0 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_Gt0 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_Gt0 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_Gt0 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_Gt0 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_Gt0 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var Gt0 
if(V_ID == 4) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_Gt1 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_Gt1 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_Gt1 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_Gt1 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_Gt1 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_Gt1 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var Gt1 
if(V_ID == 5) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_Gt2 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_Gt2 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_Gt2 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_Gt2 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_Gt2 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_Gt2 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var Gt2 
if(V_ID == 6) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_beta0 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_beta0 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_beta0 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_beta0 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_beta0 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_beta0 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_xx<pw,pencils,pencil_sz>(deriv_evars->grad2_0_0_beta0 + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_yy<pw,pencils,pencil_sz>(deriv_evars->grad2_1_1_beta0 + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_zz<pw,pencils,pencil_sz>(deriv_evars->grad2_2_2_beta0 + Z_ID * BLK_SZ, su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_0_beta0 + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad2_0_1_beta0 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_0_2_beta0 + Z_ID * BLK_SZ , su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_1_beta0 + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_1_2_beta0 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var beta0 
if(V_ID == 7) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_beta1 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_beta1 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_beta1 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_beta1 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_beta1 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_beta1 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_xx<pw,pencils,pencil_sz>(deriv_evars->grad2_0_0_beta1 + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_yy<pw,pencils,pencil_sz>(deriv_evars->grad2_1_1_beta1 + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_zz<pw,pencils,pencil_sz>(deriv_evars->grad2_2_2_beta1 + Z_ID * BLK_SZ, su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_0_beta1 + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad2_0_1_beta1 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_0_2_beta1 + Z_ID * BLK_SZ , su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_1_beta1 + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_1_2_beta1 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var beta1 
if(V_ID == 8) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_beta2 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_beta2 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_beta2 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_beta2 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_beta2 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_beta2 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_xx<pw,pencils,pencil_sz>(deriv_evars->grad2_0_0_beta2 + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_yy<pw,pencils,pencil_sz>(deriv_evars->grad2_1_1_beta2 + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_zz<pw,pencils,pencil_sz>(deriv_evars->grad2_2_2_beta2 + Z_ID * BLK_SZ, su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_0_beta2 + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad2_0_1_beta2 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_0_2_beta2 + Z_ID * BLK_SZ , su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_1_beta2 + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_1_2_beta2 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var beta2 
if(V_ID == 9) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_B0 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_B0 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_B0 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_B0 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_B0 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_B0 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var B0 
if(V_ID == 10) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_B1 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_B1 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_B1 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_B1 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_B1 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_B1 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var B1 
if(V_ID == 11) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_B2 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_B2 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_B2 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_B2 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_B2 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_B2 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var B2 
if(V_ID == 12) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_gt0 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_gt0 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_gt0 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_gt0 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_gt0 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_gt0 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_xx<pw,pencils,pencil_sz>(deriv_evars->grad2_0_0_gt0 + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_yy<pw,pencils,pencil_sz>(deriv_evars->grad2_1_1_gt0 + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_zz<pw,pencils,pencil_sz>(deriv_evars->grad2_2_2_gt0 + Z_ID * BLK_SZ, su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_0_gt0 + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad2_0_1_gt0 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_0_2_gt0 + Z_ID * BLK_SZ , su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_1_gt0 + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_1_2_gt0 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var gt0 
if(V_ID == 13) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_gt1 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_gt1 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_gt1 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_gt1 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_gt1 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_gt1 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_xx<pw,pencils,pencil_sz>(deriv_evars->grad2_0_0_gt1 + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_yy<pw,pencils,pencil_sz>(deriv_evars->grad2_1_1_gt1 + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_zz<pw,pencils,pencil_sz>(deriv_evars->grad2_2_2_gt1 + Z_ID * BLK_SZ, su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_0_gt1 + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad2_0_1_gt1 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_0_2_gt1 + Z_ID * BLK_SZ , su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_1_gt1 + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_1_2_gt1 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var gt1 
if(V_ID == 14) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_gt2 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_gt2 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_gt2 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_gt2 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_gt2 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_gt2 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_xx<pw,pencils,pencil_sz>(deriv_evars->grad2_0_0_gt2 + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_yy<pw,pencils,pencil_sz>(deriv_evars->grad2_1_1_gt2 + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_zz<pw,pencils,pencil_sz>(deriv_evars->grad2_2_2_gt2 + Z_ID * BLK_SZ, su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_0_gt2 + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad2_0_1_gt2 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_0_2_gt2 + Z_ID * BLK_SZ , su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_1_gt2 + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_1_2_gt2 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var gt2 
if(V_ID == 15) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_gt3 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_gt3 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_gt3 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_gt3 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_gt3 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_gt3 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_xx<pw,pencils,pencil_sz>(deriv_evars->grad2_0_0_gt3 + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_yy<pw,pencils,pencil_sz>(deriv_evars->grad2_1_1_gt3 + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_zz<pw,pencils,pencil_sz>(deriv_evars->grad2_2_2_gt3 + Z_ID * BLK_SZ, su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_0_gt3 + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad2_0_1_gt3 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_0_2_gt3 + Z_ID * BLK_SZ , su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_1_gt3 + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_1_2_gt3 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var gt3 
if(V_ID == 16) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_gt4 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_gt4 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_gt4 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_gt4 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_gt4 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_gt4 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_xx<pw,pencils,pencil_sz>(deriv_evars->grad2_0_0_gt4 + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_yy<pw,pencils,pencil_sz>(deriv_evars->grad2_1_1_gt4 + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_zz<pw,pencils,pencil_sz>(deriv_evars->grad2_2_2_gt4 + Z_ID * BLK_SZ, su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_0_gt4 + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad2_0_1_gt4 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_0_2_gt4 + Z_ID * BLK_SZ , su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_1_gt4 + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_1_2_gt4 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var gt4 
if(V_ID == 17) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_gt5 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_gt5 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_gt5 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_gt5 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_gt5 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_gt5 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_xx<pw,pencils,pencil_sz>(deriv_evars->grad2_0_0_gt5 + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_yy<pw,pencils,pencil_sz>(deriv_evars->grad2_1_1_gt5 + Z_ID * BLK_SZ, su, blk);
device::__blk_deriv644_zz<pw,pencils,pencil_sz>(deriv_evars->grad2_2_2_gt5 + Z_ID * BLK_SZ, su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_0_gt5 + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad2_0_1_gt5 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_0_2_gt5 + Z_ID * BLK_SZ , su, blk);
__syncthreads();
device::__load_blk_var__<DEVICE_REAL,pw,nx>(su , deriv_evars->grad_1_gt5 + Z_ID * BLK_SZ, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad2_1_2_gt5 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var gt5 
if(V_ID == 18) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_At0 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_At0 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_At0 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_At0 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_At0 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_At0 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var At0 
if(V_ID == 19) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_At1 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_At1 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_At1 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_At1 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_At1 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_At1 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var At1 
if(V_ID == 20) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_At2 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_At2 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_At2 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_At2 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_At2 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_At2 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var At2 
if(V_ID == 21) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_At3 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_At3 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_At3 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_At3 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_At3 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_At3 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var At3 
if(V_ID == 22) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_At4 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_At4 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_At4 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_At4 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_At4 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_At4 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var At4 
if(V_ID == 23) { 
device::__load_blk_var__<DEVICE_REAL,pw,nx>   (su , var_in, blk);
device::__blk_deriv644_x<pw,pencils,pencil_sz>(deriv_evars->grad_0_At5 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_x<pw,pencils,pencil_sz>(deriv_evars->kograd_0_At5 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_y<pw,pencils,pencil_sz>(deriv_evars->grad_1_At5 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_y<pw,pencils,pencil_sz>(deriv_evars->kograd_1_At5 + Z_ID * BLK_SZ , su, blk);
device::__blk_deriv644_z<pw,pencils,pencil_sz>(deriv_evars->grad_2_At5 + Z_ID * BLK_SZ , su, blk);
device::__blk_ko_deriv42_z<pw,pencils,pencil_sz>(deriv_evars->kograd_2_At5 + Z_ID * BLK_SZ , su, blk);
return;
} // end of var At5 
