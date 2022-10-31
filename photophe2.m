function [sos1,eos1,bsos1,beos1,seaos1,bseaos1,feng]=photophe2(tran11,gpp)
%%%%%%%input:DOY of changepoints; daily GPP data
%%%%%%%ouput:time of before and after the dormant season; baseline of
%%%%%%%before and after the dormant season; start and end time, and baselinesof each growing season; number of peaks

a1=gpp;
gppmax=max(gpp);
tran1(1,1)=1;
tran1(2:length(tran11)+1,1)=tran11;
tran1(end+1,1)=365;
num=length(tran1)-1;
for i=1:num%length(tran1)-1
    ix=(tran1(i,1):1:tran1(i+1))';
    iy=a1(tran1(i,1):tran1(i+1),1);
    me1=mean(iy);
    pt=polyfit(ix,iy,1);
    k(i,1)=pt(1,1);%%slope
    me(i,1)=me1;%%%%%%mean value
    m1=find(iy==gppmax);
    if isempty(m1)==0
       mos1=[tran1(i,1),tran1(i+1)]; 
    end
    clear ix iy
end
%%%%%Identify single or multiple peaks (先识别单峰或多峰)
num1=find(diff(sign(diff(me)))==2)+1;
if me(1,1)<me(2,1) && me(end,1)<me(end-1,1)
    num2=zeros(length(num1)+2,1);
    num2(1,1)=1;
    num2(2:end-1,1)=num1;
    num2(end,1)=num;
elseif me(1,1)<me(2,1)
    num2=zeros(length(num1)+1,1);
    num2(1,1)=1;
    num2(2:end,1)=num1;
elseif me(end,1)<me(end-1,1)
    num2=zeros(length(num1)+1,1);
    num2(1:end-1,1)=num1;
    num2(end,1)=num;
else
    num2=num1;
end
s1=me(num2);
tran12(:,1)=tran1(num2,1);
tran12(:,2)=tran1(num2+1,1);
std1=std(s1);

%%%find the possible bottoms (寻找可能的波谷)
bb=(max(me)-min(me))*0.3+min(me);
b1=find(s1<bb);
b11=length(b1);
if b11<2
    tran13=[tran12(1,:);tran12(end,:)];
    s2=[s1(1,1);s1(end,1)];
else
    tran13=tran12(b1,:);
    s2=s1(b1,1);
end

num3=length(s2);
feng=num3-1;
seaos1=tran13;
bseaos1=s2;
if num3<=2%%%%%%single season
    sos1=tran13(1,:);
    eos1=tran13(2,:);
    bsos1=s2(1,:);
    beos1=s2(2,:);
else
    sos1=tran13(1:2:end,:);
    eos1=tran13(2:2:end,:);
    bsos1=s2(1:2:end,:);
    beos1=s2(2:2:end,:);
end

