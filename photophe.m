function [feng,sos1,eos1]=photophe(gpp3)
%%%%%%input: gpp daily data
%%%%%%output: number of seasons, SOS(DOY), EOS(DOY)
%%%%%%%%Repeat smoothing spline, removing outliers(重复平滑样条拟合，去除离散点)
for iii1=1:20
    %             iii1=iii1+1;
    gpp2=smoothdata(gpp3,'sgolay');
    kk=(gpp2+0.01)./(gpp3+0.01);
    san=isoutlier(kk,'grubbs');
    kk1=find(san==1);
    gpp3(kk1,1)=gpp2(kk1,1);
end
mean1=mean(kk);
std1=std(kk);
kk2=find(kk<mean1-std1);
gpp2(kk2,1)=nan;
gpp2=smoothdata(gpp2,'sgolay');
trang=findchangepts(gpp2,'MinThreshold',0.5,'MinDistance',14,'Statistic','mean');
if length(trang)<2
    trang=findchangepts(gpp2,'MaxNumChanges',2,'MinDistance',14,'Statistic','mean');
end
[sosg2,eosg2,bsos1,beos1,seaos1,bseaos1,feng]=photophe2(trang,gpp2);%%%%寻找baseline
sos1=zeros(10,3);
eos1=zeros(10,3);
feng11=feng;
if feng>1%%%%%%%%%%%%Remove false growing seasons (去除伪生长季)
    for iii=1:feng
        sosg1=seaos1(iii,1);
        eosg1=seaos1(iii+1,2);
        sea=[sosg1:1:eosg1]';
        gppsea=gpp2(sea,1);
        p1=findpeaks(gppsea);
        gppmax=max(p1);
        gppmax1(iii,1)=gppmax;
    end
    gppmax2=max(gppmax1);
    tt=0;
    for iii=1:feng
        if gppmax1(iii,1)>gppmax2*0.33
            tt=tt+1;
            num1(tt,1)=iii;
        end
    end
    if tt==1
        feng=tt;
        tt1=num1(1,1);
        seaos11=[seaos1(tt1,:);seaos1(tt1+1,:)];
        bseaos11=[bseaos1(tt1,1);bseaos1(tt1+1,1)];
    else
        seaos11=seaos1;
        bseaos11=bseaos1;
    end
else
    tt=0;
    seaos11=seaos1;
    bseaos11=bseaos1;
end
%%%%%%%%%%%Calculate the time of each growing season, including single and multiple growing seasons(计算每个生长季时间，包括单，多生长季)
for iii=1:feng
    bsos2=bseaos11(iii,1);
    beos2=bseaos11(iii+1,1);
    if feng11>1 || tt>0
        sosg1=seaos11(iii,2);
        eosg1=seaos11(iii+1,2);
    else
        sosg1=seaos11(iii,1);
        eosg1=seaos11(iii+1,2);
    end
    sea=[sosg1:1:eosg1]';
    gppsea=gpp2(sea,1);
    %%%%%%Caculate outliers
    p1=findpeaks(gppsea);
    if isempty(p1)==1
        p1=max(gppsea);
    end
    gppmax=max(p1);
    k3=find(gppsea==gppmax);
    dmax=sea(k3,1);
    %%%%%%%%gpp start time
    sth1=(gppmax-bsos2)*0.10+bsos2;%10% threshold
    sth2=(gppmax-bsos2)*0.25+bsos2;%25% threshold
    sth3=(gppmax-bsos2)*0.50+bsos2;%50% threshold
    spr=[sosg1,dmax];
    gppspr(:,1)=[sosg1:1:dmax]';
    gppspr(:,2)=gpp2(spr(1,1):spr(1,2),1);
    k31=find((gppspr(:,2)-sth1)>0);
    k311=find((gppspr(:,2)-sth1)<0);
    k32=find((gppspr(:,2)-sth2)>0);
    k321=find((gppspr(:,2)-sth2)<0);
    k33=find((gppspr(:,2)-sth3)>0);
    k331=find((gppspr(:,2)-sth3)<0);
    sos1(iii,1)=gppspr(k311(end,1),1)+1;%gppspr(k31(1,1),1);
    sos1(iii,2)=gppspr(k321(end,1),1)+1;
    sos1(iii,3)=gppspr(k331(end,1),1)+1;
    %%%%%%%%gpp end time
    sth1=(gppmax-beos2)*0.10+beos2;%10% threshold
    sth2=(gppmax-beos2)*0.25+beos2;%25% threshold
    sth3=(gppmax-beos2)*0.50+beos2;%50% threshold
    aut=[dmax,eosg1];
    gppaut(:,1)=[dmax:1:eosg1]';
    gppaut(:,2)=gpp2(aut(1,1):aut(1,2),1);
    k41=find((gppaut(:,2)-sth1)<0);
    k42=find((gppaut(:,2)-sth2)<0);
    k43=find((gppaut(:,2)-sth3)<0);
    eos1(iii,1)=gppaut(k41(1,1),1);
    eos1(iii,2)=gppaut(k42(1,1),1);
    eos1(iii,3)=gppaut(k43(1,1),1);
    
    clear gppsea gppspr gppaut k31 k32 k33 k41 k42 k43 x
end
end