classdef FilterGPUdd
    properties
        dt0
        dt1
        dt2
        %dt3
        dtPad
        d
        dt0h
        dt1h
        dt2h
        %dt3h
        dtPadh
        dh
        C0
        C1
        C3
        C4
        nBatch
        varCh
        LenRec
        maxBatch
        templateBank
        templateBank0G
        templateBank1G
        templateMaskG
        Twidth
        Ttail
        Xamp
        Nwidth
        Ntail
        Nphase
        Namp
        BiasTail
        BiasWidth
        tMx
        IndMap0
        IndMap0f
        IndMap0r
        IndMap1
        NInd0
        Nstep
        NNd
        NNstep
    end
    methods
        function obj = FilterGPUdd(nBatch,varCh,LenRec,Sampling,dt0h,dt1h,dt2h,dtPadh)
            if nargin > 4
                obj.dt0h = dt0h;
                obj.dt1h = dt1h;
                obj.dt2h = dt2h;
                %obj.dt3h = dt1h+dt0h;
                obj.dtPadh=dtPadh;
                obj.dh=2*dt0h+dt1h+dt2h+dtPadh;
                obj.dt0 = 2*obj.dt0h;
                obj.dt1 = 2*obj.dt1h;
                obj.dt2 = 2*obj.dt2h;
                obj.dtPad=2*obj.dtPadh-1;%asym.
                obj.d=2*obj.dh-1;%asym.
                obj.nBatch=nBatch;
                obj.varCh=varCh;
                obj.LenRec=floor(LenRec/5)*5;
                obj.maxBatch=floor((obj.LenRec-obj.dh)/nBatch);
                %obj.templateBankG=templateBankG;%possible spike shapes
                %obj.templateMaskG=templateMaskG;%excluded spike shapes
            else
                obj.dt0h = 300;
                obj.dt1h = 90;
                obj.dt2h = 20;
                obj.dtPadh=32;
                obj.dh=2*obj.dt0h+obj.dt1h+obj.dt2h+obj.dtPadh;
                obj.dt0 = 2*obj.dt0h;
                obj.dt1 = 2*obj.dt1h;
                obj.dt2 = 2*obj.dt2h;
                obj.dtPad=2*obj.dtPadh-1;%asym.
                obj.d=2*obj.dh-1;%asym.
                obj.nBatch=nBatch;
                obj.varCh=varCh;
                obj.LenRec=floor(LenRec/5)*5;
                obj.maxBatch=floor((obj.LenRec-obj.dh)/nBatch);
            end
            obj.Nstep=obj.nBatch/obj.dt2+1;
            obj.NNd=5;
            obj.NNstep=(obj.Nstep-1)/obj.NNd;
            obj.C0=gpuArray(hanning(obj.dt0+1)/(obj.dt0h+1));
            obj.C1=gpuArray(0.1*hanning(obj.dt1+obj.dt0+1)/(obj.dt1h+obj.dt0h+1));
            obj.C1(obj.dt0h+1:end-obj.dt0h,1)=obj.C1(obj.dt0h+1:end-obj.dt0h,1)+gpuArray(0.9*hanning(obj.dt1+1)/(obj.dt1h+1));
            obj.C3=gpuArray(hanning(91));
            obj.C4=gpuArray(hanning(181));
            %obj.C4(31:end-30)=obj.C4(31:end-30)+hanning(91);
            %obj.C3=gpuArray(ones(obj.dt3+1,1)/(obj.dt3+1));
            %last index
            obj.Twidth=0.12*2.^((0:43)/12)-0.06;%template width in ms; 30 bins, want to divide in 3 (5) blocks
            %first index
            %obj.Ttail=1-(1./(2.^linspace(-1.2,3,22)+1));%proportion of the tail
            %obj.Ttail=0.5+0.03*(-6:10);%proportion of the tail
            obj.Ttail=0.5+0.04*(-7:11);%+0.002*(-6:10).*abs(-6:10);%proportion of the tail
            obj.BiasTail=sqrt(1-0.05*(abs(1-2*obj.Ttail)./(1+abs(1-2*obj.Ttail))).^2);%-0.25*
            %obj.BiasTail=sqrt(1-(abs(1-2*obj.Ttail)./(1+abs(1-2*obj.Ttail))).^2);%larger penalty!
            %obj.BiasWidth=sqrt(1-0.1*(max(obj.Twidth-0.1,0)./obj.Twidth));%assume some amount of global correlations
            obj.BiasTail=ones(1,length(obj.Ttail));
            obj.BiasWidth=ones(1,length(obj.Twidth));
            %second index
            obj.tMx=29;%%%need to check!!!
            %Fphase=(0:5);%-42;%frames*6, need 3 frames (and ignore if max outside current frame, can do convolution)
            obj.Nphase=6;
            %peak correction
            %Tshift=(obj.Ttail-1).*(3./(16*(obj.Ttail-1).^2)-1./(pi^2*(obj.Ttail-1).^2))...
            %    +(obj.Ttail).*(3./(16*(obj.Ttail).^2)-1./(pi^2*(obj.Ttail).^2));
            Fwidth=obj.Twidth*Sampling;
            
            obj.Nwidth=length(obj.Twidth);
            obj.Ntail=length(obj.Ttail);
            obj.Namp=64;
            obj.Xamp=logspace(0,2,obj.Namp+1)*3;
            %for the hierarchical step:
            IndWidth=3:2:41;%20%excluding extreme values (large/small width)
            IndTail=3:2:17;%8
            IndPhase=1:2:5;
            Ind0Width=[1 2:2:38];%20%excluding extreme values (large/small width)
            Ind0Tail=[1 2:2:12 13];%8
            obj.IndMap0=zeros(length(IndWidth)*length(IndTail)*length(IndPhase),1,'double','gpuArray');%one-based
            obj.IndMap0f=zeros(length(IndWidth)*length(IndTail),1,'double','gpuArray');%forward,one-based
            obj.IndMap0r=zeros(length(IndWidth)*length(IndTail)*length(IndPhase),1,'double','gpuArray');%reverse
            obj.NInd0=length(IndWidth)*length(IndTail);
            obj.IndMap1=zeros(7*7*obj.Nphase,1,'double','gpuArray');%one-based
            for i=1:length(IndTail)
                for j=1:length(IndWidth)
                    obj.IndMap0f(i+(j-1)*length(IndTail),1)=...
                            (Ind0Width(j)-1)*obj.Ntail*obj.Nphase+(Ind0Tail(i)-1);
                    for k=1:length(IndPhase)
                        obj.IndMap0(i+((k-1)+(j-1)*length(IndPhase))*length(IndTail),1)=...
                            (IndWidth(j)-1)*obj.Ntail*obj.Nphase+IndTail(i)+(IndPhase(k)-1)*obj.Ntail;
                        obj.IndMap0r(i+((k-1)+(j-1)*length(IndPhase))*length(IndTail),1)=...
                            i+(j-1)*length(IndTail);
                    end
                end
            end
            %obj.IndMap0(length(IndWidth)*length(IndTail)*length(IndPhase)+1:end,1)=...
            %    obj.Nwidth*obj.Ntail*obj.Nphase+(1:16)';
            for i=1:7
                for j=1:7
                    for k=1:obj.Nphase
                        obj.IndMap1(i+((k-1)+(j-1)*obj.Nphase)*7,1)=...
                            i+(k-1)*obj.Ntail+(j-1)*obj.Ntail*obj.Nphase;
                    end
                end
            end
            N1=3000;
            %make a better cosine kernel, use interpolation
            x=[(0:N1/3-1)*pi/(N1/3) cos(pi*(-N1/6:N1/6)/(N1/3))+pi (N1/3-1:-1:0)*pi/(N1/3)];
            xw=N1/sum(x>max(x)/2);
            %first, fix the tail
            obj.templateBank=zeros(2*obj.dtPadh,obj.Nwidth*obj.Ntail*obj.Nphase);
            for i=1:obj.Ntail
                y=[linspace(0,N1/2-2*obj.Ttail(i),(2-2*obj.Ttail(i))*N1/2+0.1) ...
                    linspace(N1/2,N1,(2*obj.Ttail(i))*N1/2+1.1)];
                %warped kernel
                z=interp1(0:N1,x,y);
                %relative peak location
                zPk=sum(linspace(0,1,N1+1).*z.^2)./sum(z.^2);
                for j=1:obj.Nwidth
                    template=interp1(linspace(0,1,N1+1)-zPk,z,((-obj.Nphase+1:63*obj.Nphase+2)-obj.tMx*obj.Nphase)/...
                        (xw*Fwidth(j)*obj.Nphase),'linear',0);
                    for k=1:obj.Nphase
                        obj.templateBank(:,i+(k-1)*obj.Ntail+(j-1)*obj.Ntail*obj.Nphase)=-obj.BiasTail(i)*...
                            obj.BiasWidth(j)*template(obj.Nphase-k+1:obj.Nphase:end-k)/...
                            sqrt(sum(abs(template(obj.Nphase-k+1:obj.Nphase:end-k)).^2));%unit length
                    end
                end
            end
            %want squared max at obj.tMx=28 + Fphase
            %Erange=-obj.Nphase*obj.dtPadh+1-Fphase(end)-NX:obj.Nphase*obj.dtPadh-Fphase(1)+NX;
            %Erange=-obj.Nphase*(obj.tMx+1)+2:obj.Nphase*(2*obj.dtPadh-obj.tMx);
            %obj.templateBank=zeros(2*obj.dtPadh,obj.Nwidth*obj.Ntail*obj.Nphase+16);
            %for i=1:obj.Ntail
            %    for j=1:obj.Nwidth
            %        Xrange=Erange-obj.Nphase*(Tshift(i)*Fwidth(j));
            %        fPhase0=min(max(Xrange*pi/12/Fwidth(j)/(1-obj.Ttail(i)),-pi),0);
            %        fPhase1=max(min(Xrange*pi/12/Fwidth(j)/obj.Ttail(i),pi),0);
            %        template=(-diff(sin(fPhase0))-diff(fPhase0))*(1-obj.Ttail(i))...
            %            +(-diff(sin(fPhase1))-diff(fPhase1))*obj.Ttail(i);
            %        for k=1:obj.Nphase
            %            obj.templateBank(:,i+(k-1)*obj.Ntail+(j-1)*obj.Ntail*obj.Nphase)=obj.BiasTail(i)*...
            %                template(obj.Nphase-k+1:obj.Nphase:end-k)/...
            %                sqrt(sum(abs(template(obj.Nphase-k+1:obj.Nphase:end-k)).^2));%unit length
            %        end
            %    end
            %end
            %for i=1:8
            %    obj.templateBank(1:obj.tMx+i+8,end-16+i)=-1./sqrt(obj.tMx+i);
            %    obj.templateBank(obj.tMx-i-8:end,end-8+i)=-1./sqrt(obj.dtPad-obj.tMx-i+2);
            %end
            obj.templateBank0G=gpuArray(obj.templateBank(:,obj.IndMap0));
            obj.templateBank1G=cell(length(obj.IndMap0f),1);
            for i=1:length(obj.IndMap0f)
                obj.templateBank1G{i}=gpuArray(obj.templateBank(:,obj.IndMap0f(i)+obj.IndMap1));
            end
            obj.templateMaskG=ones(obj.Ntail,obj.Nphase,obj.Nwidth,'logical','gpuArray');
            obj.templateMaskG(:,:,1:2)=false;
            obj.templateMaskG(:,:,end-1:end)=false;
            obj.templateMaskG=reshape(obj.templateMaskG,[],1);
        end
        function [LocMaxContinous,LocMaxEvents,LocMaxHist,LocMaxHist0,xbV]=FiltSection(obj,OutFileCAR,ii,i)
            if i==1
                x=zeros(obj.nBatch+obj.d,1);
                x((obj.dh+1):end)=double(h5read(OutFileCAR,'/recordings/0/data',...
                    [ii,1],[1,obj.nBatch+obj.dh-1])')/obj.varCh(ii);
                x(1:(obj.dh),1)=x(obj.d+1:-1:(obj.dh+1),1);
                xG=gpuArray(x);
            elseif i<=obj.maxBatch
                xG=gpuArray(double(h5read(OutFileCAR,'/recordings/0/data',...
                    [ii,(i-1)*obj.nBatch+1-obj.dh],[1,obj.nBatch+obj.d])')/obj.varCh(ii));
            else
                %L=mod(obj.LenRec,obj.nBatch);
                x=zeros(obj.nBatch+obj.d,1);
                x(1:end-obj.dh+1)=double(h5read(OutFileCAR,'/recordings/0/data',...
                    [ii,obj.LenRec-obj.nBatch+1-obj.dh],[1,obj.nBatch+obj.dh])')/obj.varCh(ii);
                x(end-obj.dh+2:end,1)=x(end-obj.dh:-1:end-obj.d+1,1);
                xG=gpuArray(x);
            end
            %heal bad voltage steps (need to distinguish from continous
            %voltage drifts)
            xcG=conv(diff(xG),obj.C3,'same');
            xbsG=conv(diff(xG),obj.C4,'same');
            xbsG=(diff(xG).*xcG).*max(min(2-abs(xbsG)./(abs(xcG)+1),1),0);
            xq=find(xbsG>5000);
            if ~isempty(xq)
                a=xbsG(xq);
                xqq=xq(sum((abs(xq-xq')<30).*(a-a'<0),2)==0);
                xbsG=zeros(size(xG),'double','gpuArray');
                xbsG(xqq+1)=xcG(xqq);
                xG=xG-cumsum(xbsG);
            end
            %wide kernel to remove slow drifts
            xcG=conv(xG,obj.C0,'same');
            xG=xG-xcG;
            %mask out extreme values
            %convolve with 20 ms kernel (hanning, subtract),i.e. ~60 frames
            xbsG=conv(xG,obj.C1,'same');
            xstdG=median(std(reshape(xbsG(obj.dh+1:end-obj.dh+1),obj.nBatch/50,[]),1,1));
            %determine standard deviation (for one second or batch), kill (zero) values >3stdev
            %xmaskG=((abs(xbsG)-3*xstdG)<0);
            %convolve again, same kernel, normalize by convolution of mask
            xcG=conv(xG.*((abs(xbsG)-3*xstdG)<0),obj.C1,'same')./max(conv(((abs(xbsG)-3*xstdG)<0),obj.C1,'same'),1e-6);
            %xcG=(conv(xG.*xmaskG,obj.C0,'same')+1e-3*conv(xG.*(~xmaskG),obj.C0,'same'))./...
            %    (conv(xmaskG,obj.C0,'same')+1e-3*conv(~xmaskG,obj.C0,'same'));
            %need to fix boundary!--cut after preprocessing
            xbsG=xG-xcG;
            %stdev, kill values >3stdev
            xstdG=median(std(reshape(xbsG(obj.dh+1:end-obj.dh+1),obj.nBatch/50,[]),1,1));
            %xmaskG=(abs(xbsG)<=3*xstdG).*xmaskG;
            %convolve again (10ms, multiscale), obtain baseline, subtract
            xcG=conv(xG.*((abs(xbsG)-3*xstdG)<0),obj.C1,'same')./max(conv(((abs(xbsG)-3*xstdG)<0),obj.C1,'same'),1e-6);
            %xcG=(conv(xG.*xmaskG,obj.C1,'same')+1e-3*conv(xG.*(~xmaskG),obj.C1,'same'))./...
            %    (conv(xmaskG,obj.C1,'same')+1e-3*conv(~xmaskG,obj.C1,'same'));
            %xcG=cumsum(xG(dt0+1:end-dt0,:).*xmaskG,1);
            %xcmG=cumsum(xmaskG,1);
            xbsG=xG-xcG;
            %get variance measure for intervals of 40 ms or so.
            if i>obj.maxBatch
                n=obj.LenRec-obj.maxBatch*obj.nBatch;
                xbV=gather(mean(reshape(abs(xbsG(obj.dh+obj.nBatch-n+1:obj.dh+obj.nBatch-n+...
                    obj.dt0*(ceil(n/obj.dt0)))),obj.dt0,[]),1));
                xbV=[(xbV(1:end-1)+xbV(2:end))/2 xbV(end) xbV(end)];%make one bin longer
            else
                xbV=gather(mean(reshape(abs(xbsG(obj.dh-obj.dt0h+1:end-obj.dh+obj.dt0h+1)),obj.dt0,[]),1));
                xbV=(xbV(1:end-1)+xbV(2:end))/2;
            end
            PwrG=zeros(obj.nBatch+obj.dt2,1,'double','gpuArray');
            IndG=zeros(obj.nBatch+obj.dt2,1,'double','gpuArray');
            for q=1:obj.Nstep:obj.nBatch
                P0G=zeros(obj.Nstep,1,'double','gpuArray');
                I0G=ones(obj.Nstep,1,'double','gpuArray');
                IG1=cell(obj.NNd,1);
                h=reshape(repmat(xbsG(obj.dt0+obj.dt1h+q:obj.dt0+obj.dt1h+q+obj.dtPad+obj.Nstep+1,1),1,obj.dtPad+2),[],1);
                h=reshape(h(1:(obj.dtPad+obj.Nstep+3)*(obj.dtPad+1)),(obj.dtPad+obj.Nstep+3),obj.dtPad+1);
                for qq=1:obj.NNd-1
                    [PG0,IG0]=max(h((qq-1)*obj.NNstep+1:qq*obj.NNstep,:)*obj.templateBank0G,[],2);
                    IG1{qq}=obj.IndMap0r(IG0).*(PG0>0);
                end
                [PG0,IG0]=max(h((obj.NNd-1)*obj.NNstep+1:obj.NNd*obj.NNstep+1,:)*obj.templateBank0G,[],2);
                IG1{obj.NNd}=obj.IndMap0r(IG0).*(PG0>0);
                %IG1{obj.NNd+1}=zeros(2,1);
                IG1=cat(1,IG1{:});
                for qq=1:obj.NInd0
                    XG=find(IG1==qq);
                    [PG,IG]=max(h(XG,:)*obj.templateBank1G{qq},[],2);
                    P0G(XG,1)=PG;
                    I0G(XG,1)=obj.IndMap1(IG)+obj.IndMap0f(qq);
                end
                PwrG(q:q+obj.Nstep-1,1)=P0G;
                IndG(q:q+obj.Nstep-1,1)=I0G;
            end
            %kill stuff that has a peak somewhere else
            %MaskG=obj.templateMaskG(IndG);
            %MaskG(1:obj.dtPadh)=false;
            %MaskG(end-obj.dtPadh+2:end)=false;
            %cut out relevant part? -- don't want to reshape!
            %PwrG=PwrG(obj.dtPadh+1:end-obj.dtPadh,1).*MaskG(obj.dtPadh+1:end-obj.dtPadh,1);
            %translate to shape parameter
            WidthG=floor((IndG-1)/(obj.Ntail*obj.Nphase));%zero-based
            TailG=mod(IndG-1,obj.Ntail)+1;%one-based
            PhaseG=mod(floor((IndG-1)/obj.Ntail),obj.Nphase);%zero-based
            %kill duplicates (first by amplitude, then apply mask)
            SpkG=find(PwrG(1:end-obj.dt2h,1)>0);
            MaskX=ones(size(IndG),'logical','gpuArray');
            %MaskX(~MaskG)=false;
            %AmpX=1*PwrG;
            %IndX=gpuArray(1:length(PwrG));
            for kk=1:obj.dt2h
                %Spk2=find(MaskG(SpkG+kk));
                Spk2=find(IndG(SpkG+kk)>1);
                Spk2=Spk2(abs(WidthG(SpkG(Spk2)+kk)-WidthG(SpkG(Spk2)))<7);
                %Spk2=Spk2(abs(TailG(SpkG(Spk2)+kk)-TailG(SpkG(Spk2)))<5);
                Spk2=Spk2((WidthG(SpkG(Spk2)+kk)+WidthG(SpkG(Spk2))-2*kk)>0);
                dAmpG=PwrG(SpkG(Spk2)+kk)>PwrG(SpkG(Spk2));
                MaskX(SpkG(Spk2)+kk)=MaskX(SpkG(Spk2)+kk).*dAmpG;
                MaskX(SpkG(Spk2))=MaskX(SpkG(Spk2)).*(~dAmpG);
            end
            PwrG=PwrG.*MaskX;
            %TailG=TailG.*MaskX;%dont't need
            %WidthG=WidthG.*MaskX;
            %PhaseG=PhaseG.*MaskX;
            %better:split time between spikes (locMax) evenly.
            %x=find(PwrG>0)
            %[1 round(x(1:end-1)+x(2:end)+PhaseG(...)) length(PwrG)]
            %somehow need to get at least one frame -->+diff(x)<1 of
            %something? posthoc change of some labels faster (easier)
            x=find(PwrG>obj.Xamp(1));
            [y,yi]=max(reshape(PwrG(obj.dt2h+obj.dtPadh-obj.tMx:end-obj.dt2h+obj.dtPadh-obj.tMx-1),5,(obj.nBatch)/5),[],1);
            [y1,yi1]=max(cat(1,100*y,[0 y(1:end-1)],[y(2:end) 0]),[],1);
            yi2=find(yi1==2);
            yi(yi2)=yi(yi2-1)-5;
            yi2=find(yi1==3);
            yi(yi2)=yi(yi2+1)+5;
            z=((obj.dt2h/5:(obj.nBatch+obj.dt2h)/5-1)*5+yi+obj.dtPadh-obj.tMx-1).*(y+y1>1)+(y+y1<=1);
            LocMaxC=(z'>1).*(obj.Ntail*obj.Nwidth*floor(min(log10(max(PwrG(z)/3,1))*32,64))+...
                obj.Ntail*WidthG(z)+TailG(z));%make rejected spikes zero
            x=x(abs(x-obj.nBatch/2-obj.dt2h-0.5)<(obj.nBatch/2));%no double detections!
            LocMaxE=zeros(length(x),2,'double','gpuArray');
            LocMaxE(:,1)=obj.Ntail*obj.Nwidth*floor(min(log10(max(PwrG(x)/3,1))*32,64))+...
                obj.Ntail*WidthG(x)+TailG(x);
            LocMaxE(:,2)=(i-1)*obj.nBatch+x-obj.dtPadh-obj.dt2h+obj.tMx-1+(PhaseG(x))/obj.Nphase;%time, one-based
            %create three outputs:
            %-(Nx2)-matrix of spike times (exact), power width and tail
            %-nearest local maxima, downsampled to 1/5 of sampling rate (5kHz)
            %-shape histograms
            if i>obj.maxBatch
                LocMaxE(:,2)=LocMaxE(:,2)-i*obj.nBatch+obj.LenRec;%last bin 'LenRec'.
                LocMaxE=LocMaxE(LocMaxE(:,2)>(i-1)*obj.nBatch,:);
                LocMaxContinous=[gather(LocMaxC(end-(obj.LenRec-(i-1)*obj.nBatch)/5+1:end));0];
                LocMaxEvents=gather(LocMaxE);
                LocMaxHist=gather(histcounts(LocMaxE(:,1),1:(obj.Namp*obj.Nwidth*obj.Ntail+1))');
                LocMaxHist0=gather(histcounts(LocMaxE((xbV(floor((LocMaxE(:,2)-(i-1)*obj.nBatch)/obj.dt0)+1)'<20)>0,1)...
                    ,1:(obj.Namp*obj.Nwidth*obj.Ntail+1))');
            else
                LocMaxContinous=gather(LocMaxC);
                LocMaxEvents=gather(LocMaxE);
                LocMaxHist=gather(histcounts(LocMaxE(:,1),1:(obj.Namp*obj.Nwidth*obj.Ntail+1))');
                LocMaxHist0=gather(histcounts(LocMaxE(xbV(min(max(floor((LocMaxE(:,2)-(i-1)*obj.nBatch)/obj.dt0),0)+1,...
                    obj.nBatch/obj.dt0))<20,1),1:(obj.Namp*obj.Nwidth*obj.Ntail+1))');
            end
        end
    end
end
