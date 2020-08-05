function [L,Lpwt,Spwt,Lspk]=NspikeRad(SH,SHall,Ntail,Nwidth,Namp)
%CkerOffset=[[6,4,2,0];[6,4,2,0];[6,4,2,0]];
Cker=permute([0.5;1;0.5],[1 2 3]).*permute([0.5;1;0.5],[2 1 3]).*permute([0.5;1;0.5],[2 3 1]);
Cker=Cker*1./sum(Cker,'all');
Xker=([0.5;1;0.5].*[0.5;1;0.5]')/sum([0.5;1;0.5].*[0.5;1;0.5]','all');

Ckx=reshape([zeros(6,1);hanning(27);zeros(6,1)].*hanning(13)',[],1);
Ckx=reshape(Ckx(1:end-13),38,13);
Ckx=flipud(Ckx(18:2:end-3,2:end-1));

Cka=permute(hanning(13),[1 2 3]).*permute(Ckx,[3 2 1]);
%Ckd=((-6:6).^2)'+((-6:6).^2)+permute((2*(-6:6).^2)',[2 3 1]);
%Cka(Ckd>27)=0;
Cka=Cka(2:end-1,:,:);%/sum(Cka(3:end-2,3:end-2,3:end-2),'all');%limit extent

ClipRate=0.015;%increased resolution

%create a kernel for quartile filter for watershed algo
%A=reshape((((abs(-5:Ntail+4).^2))'+((abs(-5:Nwidth+4)).^2)+permute((2*(abs(-4:4)).^2)',[2 3 1])),[],1);
%A=(([25 16 9 3 1 0 1 3 (3:Ntail+4).^2]))'+permute([16 9 4 2 1 0 1 2 (2:Nwidth+3).^2]'+...
%    max(abs(2*(-4:4)-(-5:Nwidth+4)')-1,0).^2,[3 1 2]);
%A=(([25 16 9 4 1 0 1 4 (3:Ntail+4).^2]))'+permute([16 9 4 2 1 0 1 2 (2:Nwidth+3).^2]'+...
%    max(abs(2*(-4:4)-(-5:Nwidth+4)')-1,0).^2,[3 1 2]);
%A=min(1.00101*((abs(-5:Ntail+4).^2))'+permute(1.00001*(abs(-5:Nwidth+4).^2)'+...
%    max(abs(2*(-4:4)-(-5:Nwidth+4)')-1,0).^2,[3 1 2]),25);
A=min((abs((-5:Ntail+4)+0.01).^2)'+permute((abs((-5:Nwidth+4)-0.001).^2)'+...
    max(abs(2*(-4:4)+0.0001-(-5:Nwidth+4)')-1,0).^2,[3 1 2]),25);
%%Cka=Cka.*(A(1:11,1:11,1:9)<25);
%Cka=Cka/sum(Cka,'all');
A=reshape(A,[],1);
[IU,~,IC]=unique(A);
%IH=floor(cumsum(histcounts(IC,1:max(IC)+1))/12);
%Amx=IH(end-1);
%AX=IH(IC);
Amn=20;
Amx=length(IU)-1;
AInd=find(IC<=Amx);
[ADist,ASind]=sort(IC(AInd));
ADist=max(ADist-Amn,0);
[~,~,ADist]=unique(ADist);
Amax=max(ADist);
ADind=AInd(ASind);
[~,Aoffset]=min(A);
SH0=zeros(Ntail+10,Nwidth+10,Namp+8);

HannKL1=permute(hanning(5),[2 1 3 4]).*permute(hanning(5),[2 3 1 4]).*permute(hanning(5),[2 3 4 1]);
HannKL1=HannKL1/sum(HannKL1,'all');
HannKL2=permute(hanning(9),[2 1 3 4]).*permute(hanning(9),[2 3 1 4]).*permute(hanning(9),[2 3 4 1]);
HannKL2=HannKL2/sum(HannKL2,'all');

SHx=zeros(Ntail+4,Nwidth+4,Namp+4);
SHx(3:end-2,3:end-2,6:end-2)=SH(:,:,4:end);
SHx(1:3,3:end-2,6:end-2)=repmat(SHx(3,3:end-2,6:end-2)*0.75*...
    sum(SHx(4,4:end-3,6:end-2),'all')/max(sum(SHx(3,4:end-3,6:end-2),'all'),1e-12),3,1,1);
SHx(end-2:end,3:end-2,6:end-2)=repmat(SHx(end-2,3:end-2,6:end-2)*0.75*...
    sum(SHx(end-3,4:end-3,6:end-2),'all')/max(sum(SHx(end-2,4:end-3,6:end-2),'all'),1e-12),3,1,1);
SHx(:,1:3,6:end-2)=repmat(SHx(:,4,6:end-2),1,3,1);
SHx(:,end-2:end,6:end-2)=repmat(SHx(:,end-2,6:end-2)*0.75*...
    sum(SHx(4:end-3,end-3,6:end-2),'all')/max(sum(SHx(4:end-3,end-2,6:end-2),'all'),1e-12),1,3,1);
SHx(:,:,1:6)=repmat(convn(mean(SHx(:,:,6:7),3),Xker,'same'),1,1,6);
SHavg0=convn(SHx(2:end-1,2:end-1,2:end-1),Cka,'full');
SHavg0(6:10,:,:)=SHavg0(6:10,:,:)+SHavg0(5:-1:1,:,:);
SHavg0(end-9:end-5,:,:)=SHavg0(end-9:end-5,:,:)+SHavg0(end:-1:end-4,:,:);
SHavg0(:,6:10,:)=SHavg0(:,6:10,:)+SHavg0(:,5:-1:1,:);
SHavg0(:,end-9:end-5,:)=SHavg0(:,end-9:end-5,:)+SHavg0(:,end:-1:end-4,:);
SHavg0(:,:,5:8)=SHavg0(:,:,5:8)+SHavg0(:,:,4:-1:1);
SHavg0(:,:,end-7:end-4)=SHavg0(:,:,end-7:end-4)+SHavg0(:,:,end:-1:end-3);
SHavg=SHavg0(6:end-5,6:end-5,5:end-4);
SHthr=ones(Ntail+10,Nwidth+10,Namp+8)*ClipRate;
SHthr(5:end-4,5:end-4,4:end-3)=ClipRate+0.7*(SHavg).^0.9;%%0.85
SHthr=reshape(SHthr,[],1);
SH0(4:end-3,4:end-3,3:end-2)=SHx;
SH1=reshape(SH0,[],1);
SH2c=zeros(length(SH1),1);
SH2v=zeros(length(SH1),1);
for jj=1:Amax
    SH2inc=zeros(length(SH1),1);
    Ind=find(ADist==jj);
    for j=1:length(Ind)
        SH2inc(Aoffset:end-Aoffset,1)=SH2inc(Aoffset:end-Aoffset,1)+SH1(ADind(Ind(j)):end-2*Aoffset+ADind(Ind(j)));
    end
    ThrI=find(((SH2c-SHthr)<0).*((SH2c+SH2inc-SHthr)>=0));
    if jj==1
        SH2v(ThrI,1)=jj;
    else
        SH2v(ThrI,1)=jj-(SH2c(ThrI,1)+SH2inc(ThrI,1)-SHthr(ThrI,1))./SH2inc(ThrI,1);
    end
    SH2c(Aoffset:end-Aoffset,1)=SH2c(Aoffset:end-Aoffset,1)+SH2inc(Aoffset:end-Aoffset,1);
end
ThrI=(SH2c-SHthr)<0;
SH2v(ThrI,1)=Amax;
SH2=reshape(SH2v,Ntail+10,Nwidth+10,Namp+8);
SH3=reshape(SH2(5:end-4,5:end-4,4:end-3),Ntail+2,Nwidth+2,Namp+2);
%smoothen and discretize
SH4=zeros(Ntail,Nwidth,Namp,27);
for jj=1:27
    j1=mod((jj-1),3);
    j2=mod(floor((jj-1)/3),3);
    j3=floor((jj-1)/9);
    SH4(:,:,:,jj)=SH3(j1+1:end-2+j1,j2+1:end-2+j2,j3+1:end-2+j3);
end
SH3=floor(squeeze(median(SH4,4))/2)*2;
%SH3=squeeze(median(SH4,4));
%SH3=floor(convn(SH3,Cker,'valid')*5)/5;
L0 = watershed(SH3,18);
L0(SH3==Amax)=0;
Ln=max(L0,[],'all');

%want to avoid 'narrow' clusters
Lstacked=zeros(Ln,Ntail,Nwidth,Namp);
for j=1:Ln
    Lstacked(j,:,:,:)=(L0==j);
end
Lstacked1=convn(Lstacked,HannKL1,'same')>0.8;
Lstacked2=convn(Lstacked1,HannKL2,'same');
%make sure at least 50% remaining
Lvalid=ones(Ln,1,'logical');
for j=1:Ln
    Overlap=sum(reshape(Lstacked2.*(Lstacked(j,:,:,:)>0),Ln,[]),2);
    if Overlap(j)<0.5*sum(Overlap)
        Lvalid(j,1)=false;
    end
end

[~,L]=max(Lstacked2(Lvalid,:,:,:),[],1);
L=squeeze(L.*(sum(Lstacked2,1)>0));
if ~isempty(L)
    Ln=max(L,[],'all');
else
    Ln=0;
end
%
[XW,XT,XA] = meshgrid(1:Nwidth,1:Ntail,1:Namp);
Spwt=zeros(Ln,3);
Lpwt=zeros(Ln,3);
Lspk=zeros(Ln,1);
for kk=1:Ln
    Lspk(kk,1)=sum(SHall.*(L==kk),'all');
    Lpwt(kk,1)=sum(SHall.*(L==kk).*XA,'all')/Lspk(kk,1);
    Lpwt(kk,2)=sum(SHall.*(L==kk).*XW,'all')/Lspk(kk,1);
    Lpwt(kk,3)=sum(SHall.*(L==kk).*XT,'all')/Lspk(kk,1);
    Spwt(kk,1)=sqrt(sum(SHall.*(L==kk).*XA.^2,'all')/Lspk(kk,1)-Lpwt(kk,1)^2);
    Spwt(kk,2)=sqrt(sum(SHall.*(L==kk).*XW.^2,'all')/Lspk(kk,1)-Lpwt(kk,2)^2);
    Spwt(kk,3)=sqrt(sum(SHall.*(L==kk).*XT.^2,'all')/Lspk(kk,1)-Lpwt(kk,3)^2);
end
end