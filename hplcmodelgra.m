% hplcmodelgra(Param,tg,fio,fif,to,td,te,linspace(0,tg(j)+10,2000))
function [tr] = hplcmodelgra(Param,tg,fio,fif,to,td,te,time)

fi = fio - (fio-fif)./tg.*(time-td);
fi(time<td) = fio;
fi(time>tg+td) = fif;

logkw   = squeeze(Param(:,:,1));											
logka   = squeeze(Param(:,:,2));	
logS2A  = squeeze(Param(:,:,3));
 
S1 = (logkw - logka).*(1+10.^logS2A);

tr = zeros(size(logkw));

for i=1:size(logkw,1)  
    for j=1:size(logkw,2) 

        logki = logkw(i,j) - S1(i,j) .* fi./ (1 + 10.^logS2A(i,j) .* fi);
        ki = 10.^logki; 
        inv_k_i = 1 ./ to ./ ki;
        cumtr  = cumtrapz(time,inv_k_i);
 
try
    if cumtr(end) >= 1 % if eluted before pump program ends
        trprim = interp1(cumtr,time,1);
    else                 % if eluted after pump programs ends
        trprim = (1 - cumtr(end)) .* ki(end) .* to + time(end);
    end
    tr(i,j) = trprim + to + te;
catch
    tr(i,j) = NaN;    
end

end
end
