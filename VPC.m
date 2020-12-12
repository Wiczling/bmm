function VPC(DV_VPC, DV, TT, Prior)

if ~Prior
DV_VPC(DV_VPC>2)=NaN;
DV_VPC(DV_VPC<-1)=NaN;
 
DV(DV>2)=NaN;
DV(DV<-1)=NaN;
end

DV = DV(:);
TT = TT(:);
   
binTT = [0.05 0.1:0.1:0.9];
uTT =   [0 diff(binTT)./2+binTT(1:end-1) 0.95];

hold on 

for i=1:length(uTT)-1
    
    ktore = (TT >= uTT(i) & TT<uTT(i+1));
  
    DV(ktore)=DV(ktore);
    
    prc_5(i)=prctile(DV(ktore),10,1); % data
    prc_50(i)=prctile(DV(ktore),50,1); % data
    prc_95(i)=prctile(DV(ktore),90,1); % data
    
    prc_vpc_05=prctile(DV_VPC(ktore,:),10,1); % for each simulated record
    prc_vpc_50=prctile(DV_VPC(ktore,:),50,1);
    prc_vpc_95=prctile(DV_VPC(ktore,:),90,1);   
    
    d_prc_vpc_05 = prctile(prc_vpc_05,[5,95],2);
    d_prc_vpc_50 = prctile(prc_vpc_50,[5,95],2);
    d_prc_vpc_95 = prctile(prc_vpc_95,[5,95],2);
     
    if ~any(isnan([uTT(i),d_prc_vpc_95(:,1), uTT(i+1)-uTT(i),d_prc_vpc_95(:,2)]))
    rectangle('Position',[uTT(i)+0.25*(uTT(i+1)-uTT(i)),d_prc_vpc_95(:,1), 0.5*(uTT(i+1)-uTT(i)),d_prc_vpc_95(:,2)-d_prc_vpc_95(:,1)],'FaceColor',[0.8705 0.92156 0.98039],'EdgeColor',[0.8705 0.92156 0.98039])
    end
    if ~any(isnan([uTT(i),d_prc_vpc_05(:,1), uTT(i+1)-uTT(i),d_prc_vpc_05(:,2)]))
    rectangle('Position',[uTT(i)+0.25*(uTT(i+1)-uTT(i)),d_prc_vpc_05(:,1), 0.5*(uTT(i+1)-uTT(i)),d_prc_vpc_05(:,2)-d_prc_vpc_05(:,1)],'FaceColor',[0.8705 0.92156 0.98039],'EdgeColor',[0.8705 0.92156 0.98039])
    end  
    
    if ~any(isnan([uTT(i),d_prc_vpc_50(:,1), uTT(i+1)-uTT(i),d_prc_vpc_50(:,2)]))
    rectangle('Position',[uTT(i),d_prc_vpc_50(:,1),(uTT(i+1)-uTT(i)),d_prc_vpc_50(:,2)-d_prc_vpc_50(:,1)],'FaceColor',[0.803 0.878 0.968],'EdgeColor',[0.803 0.878 0.968])    
    end

    if (d_prc_vpc_50(:,2)>d_prc_vpc_95(:,1))
     rectangle('Position',[uTT(i)+0.25*(uTT(i+1)-uTT(i)),d_prc_vpc_95(:,1), 0.5*(uTT(i+1)-uTT(i)),d_prc_vpc_50(:,2)-d_prc_vpc_95(:,1)],'FaceColor',[0.729 0.831 0.956],'EdgeColor',[0.729 0.831 0.956])
    end    
    if d_prc_vpc_50(:,1)<d_prc_vpc_05(:,2) 
     rectangle('Position',[uTT(i)+0.25*(uTT(i+1)-uTT(i)),d_prc_vpc_50(:,1), 0.5*(uTT(i+1)-uTT(i)),d_prc_vpc_05(:,2)-d_prc_vpc_50(:,1)],'FaceColor',[0.729 0.831 0.956],'EdgeColor',[0.729 0.831 0.956])
    end
    if d_prc_vpc_95(:,1)<d_prc_vpc_05(:,2)
      rectangle('Position',[uTT(i)+0.25*(uTT(i+1)-uTT(i)),d_prc_vpc_95(:,1), 0.5*(uTT(i+1)-uTT(i)),-d_prc_vpc_95(:,1)+d_prc_vpc_05(:,2)],'FaceColor',[0.729 0.831 0.956]-0.1,'EdgeColor',[0.729 0.831 0.956]-0.1)
    end
    
end
plot(TT,DV,'.','Color',[0.5 0.5 0.5],'MarkerSize',5)

plot(binTT,prc_5,'--ok','LineWidth',2)
plot(binTT,prc_50,'-ok','LineWidth',2)
plot(binTT,prc_95,'--ok','LineWidth',2)

set(gca,'Yscale','lin')

ylabel('logk')
xlabel('$$\varphi$$ (ACN)','Interpreter','latex')

