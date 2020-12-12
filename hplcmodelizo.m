function logkfi = hplcmodelizo(fi, Param)
    
 logkw   = squeeze(Param(:,:,1));											
 logka   = squeeze(Param(:,:,2));	
 logS2A  = squeeze(Param(:,:,3));
 
 S1 = (logkw - logka).*(1+10.^logS2A);
 logkfi = logkw - S1 .* fi./ (1 + 10.^logS2A .* fi);
