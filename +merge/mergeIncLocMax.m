function mergeIncLocMax(mBase,outBase,OutBaseCAR,File,subject)
%want to save a copy of the original cluster file, add another session,
%save new cluster file
z0=load([outBase subject '_' 'All_cat.mat']);
z=z0.z;
copyfile([outBase subject '_' 'All_cat.mat'],[outBase subject '_' 'All_cat.bak']);%move?

Apenalty=hanning(2*z.Nampshift+5);
Wpenalty=hanning(2*z.Nwshift+5);
Tpenalty=hanning(2*z.Ntshift+5);
z.penalty=1-Apenalty(3:end-2,1)'.*reshape((Tpenalty(3:end-2,1)*Wpenalty(3:end-2,1)'),...
    1,1,2*z.Ntshift+1,2*z.Nwshift+1);
z.Nrec=z.Nrec+1;
Iu=z.nUnits(1,1);

ii=z.Nrec;
%append another session to the data
m0=load([mBase  subject '_' File{:} '.mat']);
%iSess=mod(ii-1,nMerge)+1;
m=m0.m;
for qy=1:length(m.Channel)
    %want to ignore some stuff around the peak
    m.zSH{qy,1}=max(m.SH{qy}-z.fracSpkOffset*sum(m.SH{qy},'all'),0);
    m.zSH{qy,1}=m.zSH{qy}/max(sum(m.zSH{qy},'all'),z.Bsc);
end
%session dependent
z.SHist{ii}=m.Hist;
z.NullHisto{ii}=m.NullHisto;
z.RecId{ii}=m.RecId;
z.RawFile{ii}=m.RawFile;
z.nUnits(ii,1)=m.nUnits;
z.Units{ii}=zeros(m.nUnits,1);
z.Sessions(:,ii)=false;
z.LenRec(ii)=m.LenRec;
for c=1:z.Nch
    NewClust=find(m.Channel==c);
    %get a fully crossed JSD matrix, include moderate amplitude shifts
    Nx=length(NewClust);
    if Nx>0
        OldClust=find((z.Channel==c).*(z.UnitAge(:,ii-1)<z.nMerge));
        Ny=length(OldClust);
        NewInd=zeros(Nx,1);
        if Ny>0
            for qx=1:Nx
                JSD=zeros(Ny,2*z.Nampshift+1,2*z.Nwshift+1,2*z.Ntshift+1);
                for qy=1:Ny
                    for da=-z.Nampshift:z.Nampshift
                        daMin1=max(da,0);
                        daMax1=min(da,0);
                        daMin2=max(-da,0);%shift second matrix in opposite direction
                        daMax2=min(-da,0);
                        for dw=-z.Nwshift:z.Nwshift
                            dwMin1=max(dw,0);
                            dwMax1=min(dw,0);
                            dwMin2=max(-dw,0);%shift second matrix in opposite direction
                            dwMax2=min(-dw,0);
                            for dt=-z.Ntshift:z.Ntshift
                                dtMin1=max(dt,0);
                                dtMax1=min(dt,0);
                                dtMin2=max(-dt,0);%shift second matrix in opposite direction
                                dtMax2=min(-dt,0);
                                JD=0.5*max(m.zSH{NewClust(qx),1}(1+dtMin1:end+dtMax1,...
                                    1+dwMin1:end+dwMax1,1+daMin1:end+daMax1)+...
                                    z.zSH{OldClust(qy),1}(1+dtMin2:end+dtMax2,...
                                    1+dwMin2:end+dwMax2,1+daMin2:end+daMax2),z.Bsc);
                                JSD0=sum(m.zSH{NewClust(qx),1}(1+dtMin1:end+dtMax1,...
                                    1+dwMin1:end+dwMax1,1+daMin1:end+daMax1).*...
                                    log(max(m.zSH{NewClust(qx),1}(1+dtMin1:end+dtMax1,...
                                    1+dwMin1:end+dwMax1,1+daMin1:end+daMax1),m.Bsc)./JD),'all');
                                JSD1=sum(z.zSH{OldClust(qy),1}(1+dtMin2:end+dtMax2,...
                                    1+dwMin2:end+dwMax2,1+daMin2:end+daMax2).*...
                                    log(max(z.zSH{OldClust(qy),1}(1+dtMin2:end+dtMax2,...
                                    1+dwMin2:end+dwMax2,1+daMin2:end+daMax2),m.Bsc)./JD),'all');
                                JSD(qy,da+z.Nampshift+1,dw+z.Nwshift+1,dt+z.Ntshift+1)=squeeze((JSD0+JSD1)/2);
                            end
                        end
                    end
                end
                hh=min(reshape(JSD+repmat(z.penalty,Ny,1,1,1)*z.Jthreshold,Ny,...
                    (2*z.Nampshift+1)*(2*z.Nwshift+1)*(2*z.Ntshift+1)),[],2);
                [hz,hi]=min(hh,[],1);
                if  hz<z.Jthreshold
                    NewInd(qx,1)=OldClust(hi);
                end
            end
            %merge clusters with existing cluster when satisfying strict criterion
            for qx=1:Nx
                qy=NewInd(qx,1);
                if qy>0
                    if ~z.Sessions(qy,ii)
                        %keep half of previous histo
                        z.HistC{qy,1}=z.HistC{qy,1}+m.Hist{NewClust(qx),1};
                        z.Hist{ii}{qy}=m.Hist{NewClust(qx),1};
                        z.SH{qy,1}=0.5*z.SH{qy,1}+m.SH{NewClust(qx),1};
                        z.Sessions(qy,ii)=true;
                        z.Spikes{ii}{qy}=m.Spikes{NewClust(qx),1};
                        z.Times{ii}{qy}=m.Times{NewClust(qx),1};
                        z.nSpk(ii,qy,1)=m.nSpk(NewClust(qx),1);
                        z.fracBorder(ii,qy)=m.fracBorder(NewClust(qx),1);
                        z.nISI(ii,qy,:)=m.nISI(NewClust(qx),:);
                        z.rISI(ii,qy,:)=m.rISI(NewClust(qx),:);
                        z.pwt(ii,qy,:)=m.pwt(NewClust(qx),:);
                        z.Dpwt(ii,qy,:)=m.Spwt(NewClust(qx),:);
                        z.nSpkSum(qy,1)=m.nSpk(NewClust(qx),1);
                        z.shapeMask{ii}{qy}=m.shapeMask{NewClust(qx),1};
                    else
                        %previous histo already halved
                        z.HistC{qy,1}=z.HistC{qy,1}+m.Hist{NewClust(qx),1};
                        z.Hist{ii}{qy}=z.Hist{ii}{qy}+m.Hist{NewClust(qx),1};
                        z.SH{qy,1}=z.SH{qy,1}+m.SH{NewClust(qx),1};
                        tempTimes=[z.Times{ii}{qy};m.Times{NewClust(qx),1}];
                        tempSpikes=[z.Spikes{ii}{qy};m.Spikes{NewClust(qx),1}];
                        [~,sortTimes]=sort(tempTimes);
                        z.Spikes{ii}{qy}=tempSpikes(sortTimes,1);
                        z.Times{ii}{qy}=tempTimes(sortTimes,1);
                        %need to average shape parameter
                        z.shapeMask{ii}{qy}=(m.shapeMask{NewClust(qx),1} | z.shapeMask{ii}{qy});
                        z.fracBorder(ii,qy)=(z.fracBorder(ii,qy).*z.nSpk(ii,qy,1)+...
                            m.fracBorder(NewClust(qx),1).*m.nSpk(NewClust(qx),1))/...
                            max(z.nSpk(ii,qy,1)+m.nSpk(NewClust(qx),1),1e-6);
                        z.nISI(ii,qy,:)=z.nISI(ii,qy,:)+permute(m.nISI(NewClust(qx),:),[3 1 2]);
                        z.rISI(ii,qy,:)=z.rISI(ii,qy,:)+permute(m.rISI(NewClust(qx),:),[3 1 2]);
                        z.pwt(ii,qy,:)=(z.pwt(ii,qy,:).*z.nSpk(ii,qy,1)+...
                            reshape(m.pwt(NewClust(qx),:).*m.nSpk(NewClust(qx),1),1,1,3))/...
                            max(z.nSpk(ii,qy,1)+m.nSpk(NewClust(qx),1),1e-6);
                        z.Dpwt(ii,qy,:)=(z.Dpwt(ii,qy,:).*z.nSpk(ii,qy,1)+...
                            +reshape(m.Spwt(NewClust(qx),:).*m.nSpk(NewClust(qx),1),1,1,3))/...
                            max(z.nSpk(ii,qy,1)+m.nSpk(NewClust(qx),1),1e-6);
                        z.nSpk(ii,qy,1)=z.nSpk(ii,qy,1)+m.nSpk(NewClust(qx),1);
                        z.nSpkSum(qy,1)=z.nSpkSum(qy,1)+m.nSpk(NewClust(qx),1);
                    end
                    z.Units{ii}(NewClust(qx),1)=qy;
                    z.Units2sessUnits{qy,1}=[z.Units2sessUnits{qy,1}; [NewClust(qx) ii]];
                end
            end
            %update old clusters
            xx=find((z.Channel==c).*z.Sessions(:,ii));
            for qy=1:length(xx)
                z.zSH{xx(qy),1}=max(z.SH{xx(qy),1}-z.fracSpkOffset*sum(z.SH{xx(qy),1},'all'),0);
                z.zSH{xx(qy),1}=z.zSH{xx(qy),1}/max(sum(z.zSH{xx(qy),1},'all'),z.Bsc);
                z.rSH{xx(qy),1}=z.SH{xx(qy),1}/max(sum(z.SH{xx(qy),1},'all'),z.Bsc);
            end
        end
        %otherwise, create new clusters
        Nnew=sum(NewInd(:,1)==0);
        xx=find(NewInd(:,1)==0);
        for qx=1:length(xx)
            Iu=Iu+1;
            z.Channel(Iu,1)=c;
            z.HistC{Iu,1}=m.Hist{NewClust(xx(qx)),1};
            z.Hist{ii}{Iu}=m.Hist{NewClust(xx(qx)),1};
            z.SH{Iu,1}=m.SH{NewClust(xx(qx)),1};
            z.zSH{Iu,1}=m.zSH{NewClust(xx(qx)),1};
            z.rSH{Iu,1}=m.rSH{NewClust(xx(qx)),1};
            z.shapeMask{ii}{Iu}=m.shapeMask{xx(qx),1};
            z.Sessions(Iu,ii)=true;
            z.Units{ii}(NewClust(xx(qx)),1)=Iu;
            z.Units2sessUnits{Iu,1}=[NewClust(xx(qx)) ii];
            z.Spikes{ii}{Iu}=m.Spikes{NewClust(xx(qx)),1};
            z.Times{ii}{Iu}=m.Times{NewClust(xx(qx)),1};
            z.nSpk(ii,Iu,1)=m.nSpk(NewClust(xx(qx)),1);
            z.fracBorder(ii,Iu)=m.fracBorder(xx(qx),1);
            z.nISI(ii,Iu,:)=m.nISI(xx(qx),:);
            z.rISI(ii,Iu,:)=m.rISI(xx(qx),:);
            z.pwt(ii,Iu,:)=m.pwt(NewClust(xx(qx)),:);
            z.Dpwt(ii,Iu,:)=m.Spwt(NewClust(xx(qx)),:);
        end
        z.UnitAge=[z.UnitAge;ones(Nnew,ii-1)*z.nMerge];
    end
end
%change age for unobserved clusters
z.UnitAge(:,ii)=z.UnitAge(:,ii-1)+1;
z.UnitAge(z.Sessions(:,ii),ii)=1;
%get spike shapes
if strcmp(subject,'Jo')
    OutFileCAR=[OutBaseCAR z.RawFile{ii}(1:end-8) 'CAR.kwd'];
else
    OutFileCAR=[OutBaseCAR z.RawFile{ii}];
end
for j=1:z.Nch
    ChClust=find((z.Channel==j).*(z.Sessions(:,ii)));
    if ~isempty(ChClust)
        Raw=double(h5read(OutFileCAR,'/recordings/0/data',[j,1],[1,z.LenRec(ii)])');
        for k=1:length(ChClust)
            t=round(z.Times{ii}{ChClust(k)});
            t=t(t>z.NcutPre+1);
            t=t(t<z.LenRec(ii)-z.Ncut);
            TempSh=zeros(length(t),z.Ncut);
            for tt=1:z.Ncut
                TempSh(:,tt)=Raw(t-z.NcutPre+tt,1);
            end
            TempSh=TempSh-repmat(median(TempSh,2),1,z.Ncut);
            z.Shapes(ii,ChClust(k),:)=median(TempSh,1);
            z.ShapesD(ii,ChClust(k),:,:)=prctile(TempSh,[25 75],1);
        end
    end
end
z.nUnits=Iu;
for i=1:z.Nrec
    if length(z.Hist{i})<z.nUnits
        z.Hist{i}{z.nUnits}=[];
        z.Spikes{i}{z.nUnits}=[];
        z.Times{i}{z.nUnits}=[];
    end
end
save([outBase subject '_' 'All_cat.mat'],'z','-v7.3')
end