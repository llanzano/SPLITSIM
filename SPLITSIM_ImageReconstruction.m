
function []=SPLITSIM_ImageReconstruction(filename );

%                                                                   
%  Luca Lanzanò & Diaspro-Lab                  
%  Istituto Italiano di Tecnologia - Nanoscopy Department      
%  University of Catania - Department of Physics and Astronomy    
%
%  User-Friendly Version based on the article:                     
%  Chromatin investigation in the nucleus using a phasor approach to structured illumination microscopy 
%
%  by: Isotta Cainero, Elena Cerutti, Mario Faretta, Gaetano Ivan Dellino, Pier Giuseppe Pelicci, Paolo Bianchini, 
%     Giuseppe Vicidomini, Alberto Diaspro and Luca Lanzanò
%  Biopysical Journal
                    

%open image stack
[filename,pathname, filterindex] = uigetfile({'*.tif'},'Please select an image stack');
filenamefull = [pathname, filename];   
A=simpleSPLIT_readfiletif(filenamefull);  

% line for fairSIM data
% slicez=4;
% A=A(:,:,cat(2,5*(slicez-1)+1:5*(slicez-1)+5,35+(5*(slicez-1)+1:5*(slicez-1)+5),35+35+( 5*(slicez-1)+1:5*(slicez-1)+5)));

X0=size(A,1);
Y0=size(A,2);
ROI0=min(X0,Y0);
%parameters
prompt = {'Number of phase steps','Number of orientations','Oversampling','ROI'};
dlg_title = 'Input parameters'; 
num_lines = 1;
def = {'5','3','1',num2str(ROI0)};
answer = inputdlg(prompt,dlg_title,num_lines,def);

Nph=str2num(answer{1});
Nangles=str2num(answer{2});
over=str2num(answer{3});
ROI=str2num(answer{4});

Bordx=round(0.5*(size(A,1)-ROI));
Bordy=round(0.5*(size(A,2)-ROI));
A=A(1+Bordx:X0-Bordx,1+Bordy:Y0-Bordy,:);
% 
% over=1;
T=size(A,3)/Nangles;
X=size(A,1);
Y=size(A,2);

nc=2;
h=1;
%calc widefield image = sum of stack
Ntot=sum(A,3);
maxc=round(max(max(Ntot)));
max1=max(max(A(:,:,1)));

MaxNtot=max(max(Ntot));

%threshold
thr=0;
thrpc=0;
figcheck=1;
while figcheck==1 
thr=max(thr, thrpc*MaxNtot);    
B=simpleSPLIT_Threshold(Ntot,thr);
prompt2 = {'Threshold value:', 'With respect to max (0-1):'}; 
dlg_title2 = 'Input parameters'; 
num_lines = 1;
def2 = {num2str(thr),num2str(thrpc)};
answer2 = inputdlg(prompt2,dlg_title2,num_lines,def2);
figcheck=~isempty(answer2); 
if figcheck==1
thr=str2num(answer2{1});
thrpc=str2num(answer2{2});
end
end
close

%load pattern stack (for simulations only)
[filename_patt,pathname, filterindex] = uigetfile({'*.tif'},'Select the SIM pattern stack (simulations)');
if ~isequal(filename_patt,0)
flagsim=1;
max_mod=0.5;
max_mod2=0.5;
filenamefull_patt = [pathname, filename_patt];   
Ap=simpleSPLIT_readfiletif(filenamefull_patt);
Ap2=Ap;
else
flagsim=0;
T0=8.2;
T0angle=[T0 T0 T0];
alphai= [  -0.2450    0.8022   -1.2922 ];
f0=[-0.9 -0.1 0.6];
df2h=[ -0.5133   -3.9800   -1.5300 ];
max_mod=0.4236;
max_mod2=0.2754;

answer2 = questdlg('Load Pattern parameters');
if strcmp(answer2,'Yes')==1
[filenamepatt,pathnamepatt, filterindex] = uigetfile({'*_PATTERN.mat'},'Please select');
filenamefullpatt = [pathnamepatt, filenamepatt];   
load( [filenamefullpatt], 'T0', 'alphai' ,  'f0', 'df2h','max_mod', 'min_ph' , 'max_mod2', 'min_ph2' ) ;
else
prompt = {'Period','Angles (rad)','delta phase 2nd order','max mod 1', 'max mod 2' };
dlg_title = 'Pattern parameters'; 
num_lines = 1;
def = {num2str(T0),num2str(alphai),num2str(df2h), num2str(max_mod),num2str(max_mod2) };
answerp = inputdlg(prompt,dlg_title,num_lines,def);
T0=str2num(answerp{1});
alphai=str2num(answerp{2});
df2h=str2num(answerp{3});    
max_mod=str2num(answerp{4});    
max_mod2=str2num(answerp{5});    
end 

%additional phase adjustment
% ALL ANGLES harmonic 1
ntestpoints1=1;
ntestpoints2=1;
ntestpoints5=1;
d_f=0;
for ik=1
% B2=ones(X,Y);
[gf_nocorr ]=simpleSPLITsim_PhaseFromSIMPattern(A,T,1) ;
% mod_nocorr=sqrt( sum( (abs(gf_nocorr)).^2 ,3) );
mod_nocorr= abs(gf_nocorr)  ;
figcheck111=1;
thr1=thr;
thr1mod=0;
while figcheck111==1
% B1=simpleSPLIT_Threshold(Ntot,thr1).*simpleSPLIT_Threshold(mod_nocorr,thr1mod);
B1=simpleSPLIT_Threshold(Ntot,thr1).*simpleSPLIT_Threshold(mod_nocorr(:,:,1),thr1mod).*simpleSPLIT_Threshold(mod_nocorr(:,:,2),thr1mod).*simpleSPLIT_Threshold(mod_nocorr(:,:,3),thr1mod);
ecc_test=zeros(ntestpoints1,ntestpoints2,ntestpoints5);
for ntest5=1:ntestpoints5
           for ntest2=1:ntestpoints2
               for ntest1=1:ntestpoints1
                   f0_test(2)=f0(2)-d_f+(ntest5-1)*2*d_f/ntestpoints5;
                   f0_test(3)=f0(3)-d_f+(ntest2-1)*2*d_f/ntestpoints2;
                   f0_test(1)=f0(1)-d_f+(ntest1-1)*2*d_f/ntestpoints1;
                   [Ap]=SPLITsimple_simulate_pattern(100, X, Y, T0, alphai,f0_test, 1, T, 'N');
%                    [Ap1]=SPLITsimple_simulate_pattern(100, X, Y, T0, alphai(1),f0_test(1), 1, T, 'N');
%                    [Ap2]=SPLITsimple_simulate_pattern(100, X, Y, T0, alphai(2),f0_test(2), 1, T, 'N');
%                    [Ap3]=SPLITsimple_simulate_pattern(100, X, Y, T0, alphai(3),f0_test(3), 1, T, 'N');
                %calc phase reference from pattern file
                [g_ref ]=simpleSPLITsim_PhaseFromSIMPattern(Ap,T,1) ;
%                 [g_ref1 ]=simpleSPLITsim_PhaseFromSIMPattern(Ap1,T,1) ;
%                 [g_ref2 ]=simpleSPLITsim_PhaseFromSIMPattern(Ap2,T,1) ;
%                 [g_ref3 ]=simpleSPLITsim_PhaseFromSIMPattern(Ap3,T,1) ;

                gf_corr=SPLITsimple_CorrectPhasor(gf_nocorr, g_ref, B1,1) ;               
                phaseimg=(angle(conj(gf_corr))); 
                phaseimg(B1==0)=0;
                
%                 mod1_img=abs(gf_corr);
%                 mod1_img(B1==0)=0;
%                 [elem, centers]=hist(mod1_img(B1>0),10);
%                 [C,I] = max(elem);
%                 ind1 = find(elem>0.05*C, 1, 'first') ;
%                 min_mod_0=centers(ind1);
%                 ind2 = find(elem>0.05*C, 1, 'last') ;
%                 max_mod_0=centers(ind2);
%                 
%                 wlin_mod= abs(mod1_img-min_mod_0)/abs(max_mod_0-min_mod_0);  %standard
% %               wlin_mod= abs(mod1_img-0)/abs(max_mod-0);  % test
%                              
% %               wlin_mod=(1/0.13)* exp(-1./( abs(mod1_img-min_mod)./abs(max_mod-min_mod) )) ;
%                 wlin_mod=wlin_mod/max(wlin_mod, [],'all');
%                 
%                 
%                 phaseimg_w=wlin_mod.*phaseimg;  %standard
% %                 phaseimg_w=wlin_mod.*(1./phaseimg);  % test
              phaseimg_w=phaseimg;                 

                ecc=mean(phaseimg_w(phaseimg_w>0)); 
%                 wlin_mod(wlin_mod<min_mod)=0;
%                 wlin_mod(wlin_mod>max_mod)=1;
%                 gf_corr1=SPLITsimple_CorrectPhasor(gf_nocorr1, g_ref1, B) ;               
%                 phaseimg1=(angle(conj(gf_corr1))); 
%                 phaseimg1(B==0)=0;
%                 gf_corr2=SPLITsimple_CorrectPhasor(gf_nocorr2, g_ref2, B) ;               
%                 phaseimg2=(angle(conj(gf_corr2))); 
%                 phaseimg2(B==0)=0;
%                 gf_corr3=SPLITsimple_CorrectPhasor(gf_nocorr3, g_ref3, B) ;               
%                 phaseimg3=(angle(conj(gf_corr3))); 
%                 phaseimg3(B==0)=0;
%                 phaseimg=(phaseimg1+phaseimg2+phaseimg3)/3;
                
%               ecc=mean(phaseimg(phaseimg>0));                
                ecc_test(ntest1,ntest2,ntest5)=ecc;  % (angle, Phase, Period)
                end
           end
end

% find minimum and set
ecc_max=min(min(min(ecc_test)));
nindex=find(ecc_test==ecc_max,1);
siz=[ntestpoints1, ntestpoints2, ntestpoints5];
[nmax1,nmax2,nmax5] = ind2sub(siz,nindex) ;
f0(1)=f0(1)-d_f+(nmax1-1)*2*d_f/ntestpoints1;
f0(2)=f0(2)-d_f+(nmax5-1)*2*d_f/ntestpoints5; 
f0(3)=f0(3)-d_f+(nmax2-1)*2*d_f/ntestpoints2;
[Aptot]=SPLITsimple_simulate_pattern(100, X, Y, T0, alphai,f0, 1, T, 'N');
%calc phase reference from pattern file 
[g_ref ]=simpleSPLITsim_PhaseFromSIMPattern(Aptot,T,1) ;
gf_corr=SPLITsimple_CorrectPhasor(gf_nocorr, g_ref, B1,1) ; % correct for all the angles and average
phaseimg=(angle(conj(gf_corr))); 
phaseimg(B1==0)=0;
figure
% subplot(2,2,1)
% imagesc(mod1_img)
% axis image
% title(['Mod'])
% subplot(2,2,2)
imagesc(phaseimg);
axis image
title(['Phase'])
% subplot(2,2,2)
% imagesc(B.*abs(g_ns(:,:,2)));
% title(['Mod - harm=',num2str(1)])
xzoom=1+floor(X/2);
yzoom=1+floor(Y/2);
xzoom1=1+floor(X/8);
yzoom1=1+floor(Y/8);
% subplot(2,2,3)
% imagesc(phaseimg(xzoom-xzoom1:xzoom+xzoom1,yzoom-yzoom1:yzoom+yzoom1))
% axis image

prompt = {'Angle 1 phase','Angle 2 phase','Angle 3 phase','delta','test points phase','threshold','threshold_mod','Min Phase value'};
dlg_title = 'Pattern parameters'; 
num_lines = 1;
def = { num2str(f0(1)), num2str(f0(2)), num2str(f0(3)), num2str(d_f), num2str(ntestpoints5), num2str(thr1), num2str(thr1mod), num2str(ecc_max)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
figcheck111=~isempty(answer); 
if figcheck111==1
f0(1)=str2num(answer{1});
f0(2)=str2num(answer{2});
f0(3)=str2num(answer{3});
d_f=str2num(answer{4});
ntestpoints5=str2num(answer{5});
ntestpoints1=ntestpoints5;
ntestpoints2=ntestpoints5;
thr1=str2num(answer{6});
thr1mod=str2num(answer{7});
% ntestpoints1=str2num(answer{8});
% ntestpoints2=str2num(answer{9});
close all
end
end
% end of optimization of all phases
close all

end

%save new f0 value with sample filename
answer = questdlg('Save Pattern parameters?');
if strcmp(answer,'Yes')==1
    if isequal(ROI, min(X0,Y0))
save( [filename(1:end-4),'_PATTERN.mat'], 'T0', 'alphai' , 'f0', 'df2h' ,'max_mod',  'max_mod2' ) ;
    else
save( [filename(1:end-4),'_ROI',num2str(ROI),'_PATTERN.mat'], 'T0', 'alphai' , 'f0', 'df2h' ,'max_mod',  'max_mod2' ) ;
       
    end
end

%calculate final optimized pattern
[Ap]=SPLITsimple_simulate_pattern(100, X, Y, T0, alphai,f0, 1, T, 'N');
[Ap2]=SPLITsimple_simulate_pattern_h(100, X, Y, T0,2, alphai,f0+df2h, 1, T, 'N');
end

% part after optimization of pattern parameters
%calc phase reference from pattern file
for h=1:2
harm1=1+2*(h-1);

if h==1
[g_ref ]=simpleSPLITsim_PhaseFromSIMPattern(Ap,T,h) ;  % gives complex form of (g,s) for each angle at the selected harmonic
else
 [g_ref ]=simpleSPLITsim_PhaseFromSIMPattern(Ap2,T,h) ;
end

[gf_nocorr ]=simpleSPLITsim_PhaseFromSIMPattern(A,T,h) ; % idem
% gf_corr=gf_nocorr;
gf_corr(:,:,h)=SPLITsimple_CorrectPhasor(gf_nocorr, g_ref, B,h) ; % correct and average all angles
end

%calculate phasor plot (can be normalized with normf)
% normf=(1/sqrt(2))*[1/max_mod 1/max_mod2];
normf=[1 1];
sm=0;
smp=0;
[gf, HistoH]=simpleSPLITSIM_phasor_plot_corr(A(:,:,1:end),T,thr,sm,smp, gf_corr, normf) ; 
%eventualmente inserire qui le opzioni di plot con un menu
% remove normalization used in the phasor plot
for h1=1:2   
gf(:,:,h1+1)=gf(:,:,h1+1)*(1/normf(h1));
end

% gf (complex, corrected); HistoH is histo at several harmonics 
histo_ns=HistoH(:,:,1); % first harm phasor plot for export
histo_ns2=HistoH(:,:,2); % 2nd harm phasor plot for export
% without threshold
phaseimg_nothr=abs(angle(gf(:,:,2)));
mod1_nothr=abs(gf(:,:,2));
% with threshold
phaseimg=abs(angle(gf(:,:,2))).*B;
mod1=mod1_nothr;
mod1(isnan(mod1))=0;
mod1=mod1.*B;

% without threshold 2
phaseimg_nothr2=abs(angle(gf(:,:,3)));
mod2_nothr=abs(gf(:,:,3));
% with threshold
phaseimg2=abs(angle(gf(:,:,3))).*B;
mod2=mod2_nothr;
mod2(isnan(mod2) )=0;
mod2=mod2.*B;

% over=1; %oversampling if over>1
Ntot = imresize(Ntot,over);
B = imresize(B,over);
phaseimg_nothr = imresize(phaseimg_nothr,over);
mod1_nothr = imresize(mod1_nothr,over);
phaseimg_nothr2 = imresize(phaseimg_nothr2,over);
mod2_nothr = imresize(mod2_nothr,over);
phaseimg = imresize(phaseimg,over);
mod1 = imresize(mod1,over);
phaseimg2 = imresize(phaseimg2,over);
mod2 = imresize(mod2,over);

delta_ph=pi/2;
delta_ph_piunit=0.5;
min_ph=0;
min_ph2=0;
figure
figcheck3=1;
sm=0;
smr=0;
splitmode=2;
if flagsim==1
    splitmode=1;
end
while figcheck3==1 

% calculate fractions (SPLIT images) according to different formulas:
f=zeros(over*X,over*Y,2);
fbkgd=zeros(over*X,over*Y);
delta_ph2=1*delta_ph;
switch splitmode
case 1 % for sim  - only 1st harm
[f1harm,f1harmlin]=SPLITSIM_fraction_2comp(phaseimg_nothr, min_ph,delta_ph, mod1_nothr, max_mod, over);
% [f2harm,f2harmlin]=SPLITSIM_fraction_2comp(phaseimg_nothr2, min_ph2,delta_ph2, mod2_nothr, max_mod2, over);
f(:,:,1)=f1harm(:,:,1) ;
f(:,:,2)=ones(over*X,over*Y)-f(:,:,1);

case 2 % for data  - uses 2 harmonics
[f1harm,f1harmlin]=SPLITSIM_fraction_2comp(phaseimg_nothr, min_ph,delta_ph, mod1_nothr, max_mod, over);
[f2harm,f2harmlin]=SPLITSIM_fraction_2comp(phaseimg_nothr2, min_ph2,delta_ph2, mod2_nothr, max_mod2, over);
f(:,:,1)=f1harm(:,:,1).*f2harm(:,:,1) ;
f(:,:,2)=ones(over*X,over*Y)-f(:,:,1);

case 3 % for data  - only 2nd harm
[f1harm,f1harmlin]=SPLITSIM_fraction_2comp(phaseimg_nothr, min_ph,delta_ph, mod1_nothr, max_mod, over);
[f2harm,f2harmlin]=SPLITSIM_fraction_2comp(phaseimg_nothr2, min_ph2,delta_ph2, mod2_nothr, max_mod2, over);
f(:,:,1)=f2harm(:,:,1) ;
f(:,:,2)=ones(over*X,over*Y)-f(:,:,1);


end

% calculate images from fractions
for k=1:nc
SplitImg(:,:,k)=f(:,:,k).*Ntot;
end
%plot  figures
maxc=round(max(max(A(:,:,1))));
maxcg=round(max(max(Ntot)));
maxc1=round(max(max(SplitImg(:,:,1))));

for k=1:nc
SplitImg(:,:,k)=simpleSPLIT_smooth_simple(SplitImg(:,:,k),sm,smr);
end


subplot(2,4,1)
colormap(hot)
imagesc(Ntot, [maxcg-maxcg,maxcg])
title('Widefield');
axis image
subplot(2,4,2)
colormap(hot)
imagesc(SplitImg(:,:,1))
% imagesc(SplitImg(:,:,1), [maxc1-maxc1,maxc1])
title('SPLIT');
axis image
subplot(2,4,3)
colormap(hot)
imagesc(SplitImg(:,:,2))
% imagesc(SplitImg(:,:,2), [maxc1-maxc1,maxc1])
title('SPLIT out');
axis image
% subplot(2,4,4)
% colormap(hot)
% imagesc(SplitImg(:,:,3), [maxc1-maxc1,maxc1])
% title('SPLIT bkgd');
% axis image

xzoom=1+floor(over*X/2);
yzoom=1+floor(over*Y/2);
xzoom1=1+floor(over*X/8);
yzoom1=1+floor(over*Y/8);

subplot(2,4,5)
colormap(hot)
imagesc(Ntot(xzoom-xzoom1:xzoom+xzoom1,yzoom-yzoom1:yzoom+yzoom1), [maxcg-maxcg,maxcg])
title('Widefield');
axis image

subplot(2,4,6)
colormap(hot)
imagesc(SplitImg(xzoom-xzoom1:xzoom+xzoom1,yzoom-yzoom1:yzoom+yzoom1,1) )
title('SPLIT');
axis image
subplot(2,4,7)
colormap(hot)
imagesc(SplitImg(xzoom-xzoom1:xzoom+xzoom1,yzoom-yzoom1:yzoom+yzoom1,2) )
title('SPLIT out');
axis image
%end here SPLIT operations 

prompt3 = {'Phase max (Pi unit)','Smoothing factor','Smoothing cycles','Harmonics for SPLIT'};
dlg_title3 = 'SPLIT parameters'; 
num_lines3 = 1;
def3 = {num2str(delta_ph/pi,2), num2str(sm,2), num2str(smr,2),num2str(splitmode)  };
answer3 = inputdlg(prompt3,dlg_title3,num_lines3,def3);
figcheck3=~isempty(answer3); 
if figcheck3==1
delta_ph_piunit=str2num(answer3{1});
delta_ph=pi*delta_ph_piunit;
sm=str2num(answer3{2});
smr=str2num(answer3{3});
splitmode=str2num(answer3{4});
end
end




answer = questdlg('Save data?');
if strcmp(answer,'Yes')==1

filenameout=[filenamefull(1:end-4),'_SPLITSIM'] ;

%save data in Matlab
save([filenameout,'_splitsim.mat']);

% export scaled tiff image
A1=double(SplitImg(:,:,1));
MaxValImg=max(A1,[],'all');
MinValImg=min(A1,[],'all');
Aout=(A1-MinValImg)/(MaxValImg-MinValImg);
outputFileName = [filenameout, '_min',num2str(MinValImg,2),'max',num2str(MaxValImg),'.tiff'];
delete outputFileName ;
imwrite(Aout, outputFileName);


end


end

% required functions

function A=simpleSPLIT_readfiletif(fname)

info = imfinfo(fname);
nslice = numel(info);

A=imread(fname, 1);  % read tif stack data
for k = 2:nslice
    B = imread(fname, k);
    A=cat(3,A,B);
end


end

function y=simpleSPLIT_smooth_simple(M,sm,n)
y=M;
if sm>0
filt = (1/(8+1/sm))*[1 1 1; 1 1/sm 1; 1 1 1]; % sm factor <=1 
    for i=1:n
    y = filter2(filt,y);
    end
end
    
end

function B=simpleSPLIT_Threshold(A,thr)
  
if length(thr)==1
B=A;
B(B<=thr)=0;
B(B>thr)=1;
else
B=A;
B(B>thr(2))=0;
B(B<=thr(1))=0;
B(B>0)=1;
end

% figure
subplot(1,2,1)
imagesc(A,[min(nonzeros(A)),max(nonzeros(A))])
axis image
subplot(1,2,2)
imagesc(B)
colormap(hot)
axis image


end

function [g_ref ]=simpleSPLITsim_PhaseFromSIMPattern(Afull,T,h1)

  % phaseref(ix,j,angle) or T0, alphai,phase0
% h1=1;
% T is number of phase shifts, typically 5
Afull=double(Afull);
dimimg=size(Afull);
X=dimimg(1);
Y=dimimg(2);
Z=dimimg(3);
Nangles=round(Z/T);
Nph=T;
% Ntot=sum(Afull,3);
for anglesim=1:Nangles    % calculate g and s for each angle 

A=Afull(:,:,(anglesim-1)*Nph+1:anglesim*Nph);  % get substack for each angle
DC=sum(A,3);
Ntot=DC;
max1=max(max(A(:,:,1))); 
for k = 2:T
    DC=cat(3,DC,sum(A,3));
end
%caclulation of g and s
if T>4
gf=fft(A, [], 3)./DC;   %temporal fft of the stack
else
     for ix=1:X
        for j=1:Y
            for vk=1:T
   gfcos(ix,j,vk)=( sum( squeeze(A(ix,j,:))'.*cos(2*pi*(vk-1)*(0:T-1)/T)) )/DC(ix,j,1);  
   gfsin(ix,j,vk)=( sum( squeeze(A(ix,j,:))'.*sin(2*pi*(vk-1)*(0:T-1)/T)) )/DC(ix,j,1);           
   gf(ix,j,vk)=gfcos(ix,j,vk)-1i*gfsin(ix,j,vk);      
            end
        end
     end
end

test1=isfinite(gf);
test2=prod(test1,3);
for vk=1:T
    gf(:,:,vk)=gf(:,:,vk).*test2;
end

if T==2
gfsin(:,:,:)=0;
end


g_ref(:,:,anglesim)=(gf(:,:,h1+1)) ; %

end

end

function [gf, HistoH ]=simpleSPLITSIM_phasor_plot_corr(Afull,T,Thrg,sm,smr, gf_corr, normf)  % phaseref(ix,j,angle) or T0, alphai,phase0
% added normf as a factor to rescale modulus of phasor
hmax=size(gf_corr,3);
Afull=double(Afull);
dimimg=size(Afull);
X=dimimg(1);
Y=dimimg(2);
Z=dimimg(3);
Nangles=round(Z/T);
Nph=T;
Ntot=sum(Afull,3);
gf=zeros(X,Y,hmax);
for h1=1:hmax   
gf(:,:,h1+1)=normf(h1)*gf_corr(:,:,h1);
end

%smooth g and s
gf_s=gf;
for k=1:hmax          % smooth g and s
    gf(:,:,k) = simpleSPLIT_smooth_simple(gf(:,:,k),sm,smr);
    gf_s(:,:,k) = simpleSPLIT_smooth_simple(gf_s(:,:,k),sm,smr+1);
end

figure
for h1=1:hmax
mod1=abs(gf(:,:,h1+1));
mod1(isnan(mod1) | Ntot<Thrg)=0;
maxmod=max(max(mod1));
phase1=abs(angle(gf(:,:,h1+1)));
phase1(isnan(phase1) | Ntot<Thrg)=2;

%mod image for export
modimg=mod1;
modimg(isnan(modimg) | Ntot<Thrg)=0;

p=1;                          %phasor plot histo g,s  minimum smoothing
 for i=1:X
        for j=1:Y
            if Ntot(i,j)>Thrg
                gplot1(p)=conj(gf(i,j,1+h1));
                p=p+1;
            end
        end
 end

[H1, ~, ~] = simpleSPLIT_hist2d_new(gplot1,200,200,[-1 1],[1 -1]);
HistoH(:,:,h1)=H1;
subplot(hmax,3,1+3*(h1-1));
imagesc(H1);
title(['Phasor Plot - harm=',num2str(h1)])
axis image 
axis([49 200 1 151])

subplot(hmax,3,2+3*(h1-1));
if isfinite(maxmod)
imagesc(mod1,[maxmod-maxmod,maxmod]);
else
    imagesc(mod1);
end
title(['modulation'])
axis image 

subplot(hmax,3,3+3*(h1-1));
imagesc(phase1);
title(['phase'])
axis image
end

end

function [Hout Xbins Ybins] = simpleSPLIT_hist2d_new(D, varargin) %Xn, Yn, Xrange, Yrange)

    % PROCESS INPUT D
    if nargin < 1 %check D is specified
        error 'Input D not specified'
    end
    
    Dcomplex = false;
    if ~isreal(D) %if D is complex ...
        if isvector(D) %if D is a vector, split into real and imaginary
            D=[real(D(:)) imag(D(:))];
        else %Thrgow error
            error 'D must be either a complex vector or nx2 real array'
        end
        Dcomplex = true;
    end

    if (size(D,1)<size(D,2) && size(D,1)>1)
        D=D';
    end
    
    if size(D,2)~=2
        error('The input data matrix must have 2 rows or 2 columns');
    end
    
    % PROCESS OTHER INPUTS
    var = varargin;

    % check if DISPLAY is specified
    index = find(strcmpi(var,'display'));
    if ~isempty(index)
        display = true;
        var(index) = [];
    else
        display = false;
    end

    % process number of bins    
    Xn = 20; %default
    Xndefault = true;
    if numel(var)>=1 && ~isempty(var{1}) % Xn is specified
        if ~isscalar(var{1})
            error 'Xn must be scalar'
        elseif var{1}<1   %|| ~isinteger(var{1})
            error 'Xn must be an integer greater than or equal to 1'
        else
            var{1}~=floor(var{1});
            Xn = var{1};
            Xndefault = false;
        end
    end

    Yn = 20; %default
    Yndefault = true;
    if numel(var)>=2 && ~isempty(var{2}) % Yn is specified
        if ~isscalar(var{2})
            error 'Yn must be scalar'
        elseif var{2}<1   % || ~isinteger(var{2})
            error 'Xn must be an integer greater than or equal to 1'
        else
            var{2}~=floor(var{2});
            Yn = var{2};
            Yndefault = false;
        end
    end
    
    % process ranges
    if numel(var) < 3 || isempty(var{3}) %if XRange not specified
        Xrange=[min(D(:,1)),max(D(:,1))]; %default
    else
        if nnz(size(var{3})==[1 2]) ~= 2 %check is 1x2 array
            error 'XRange must be 1x2 array'
        end
        Xrange = var{3};
    end
    if Xrange(1)==Xrange(2) %handle case where XLO==XHI
        if Xndefault
            Xn = 1;
        else
            Xrange(1) = Xrange(1) - floor(Xn/2);
            Xrange(2) = Xrange(2) + floor((Xn-1)/2);
        end
    end
    
    if numel(var) < 4 || isempty(var{4}) %if XRange not specified
        Yrange=[min(D(:,2)),max(D(:,2))]; %default
    else
        if nnz(size(var{4})==[1 2]) ~= 2 %check is 1x2 array
            error 'YRange must be 1x2 array'
        end
        Yrange = var{4};
    end
    if Yrange(1)==Yrange(2) %handle case where YLO==YHI
        if Yndefault
            Yn = 1;
        else
            Yrange(1) = Yrange(1) - floor(Yn/2);
            Yrange(2) = Yrange(2) + floor((Yn-1)/2);
        end
    end
        
    % SET UP BINS
    Xlo = Xrange(1) ; Xhi = Xrange(2) ;
    Ylo = Yrange(1) ; Yhi = Yrange(2) ;
    if Xn == 1
        XnIs1 = true;
        Xbins = [Xlo Inf];
        Xn = 2;
    else
        XnIs1 = false;
        Xbins = linspace(Xlo,Xhi,Xn) ;
    end
    if Yn == 1
        YnIs1 = true;
        Ybins = [Ylo Inf];
        Yn = 2;
    else
        YnIs1 = false;
        Ybins = linspace(Ylo,Yhi,Yn) ;
    end
    
    Z = linspace(1, Xn+(1-1/(Yn+1)), Xn*Yn);
    
    % split data
    Dx = floor((D(:,1)-Xlo)/(Xhi-Xlo)*(Xn-1))+1;
    Dy = floor((D(:,2)-Ylo)/(Yhi-Ylo)*(Yn-1))+1;
    Dz = Dx + Dy/(Yn) ;
    
    % calculate histogram
    h = reshape(histc(Dz, Z), Yn, Xn);
    
    if nargout >=1
        Hout = h;
    end
    
    if XnIs1
        Xn = 1;
        Xbins = Xbins(1);
        h = sum(h,1);
    end
    if YnIs1
        Yn = 1;
        Ybins = Ybins(1);
        h = sum(h,2);
    end
    
    % DISPLAY IF REQUESTED
    if ~display
        return
    end
        
    [x y] = meshgrid(Xbins,Ybins);
    dispH = h;

    % handle cases when Xn or Yn
    if Xn==1
        dispH = padarray(dispH,[1 0], 'pre');
        x = [x x];
        y = [y y];
    end
    if Yn==1
        dispH = padarray(dispH, [0 1], 'pre');
        x = [x;x];
        y = [y;y];
    end

    surf(x,y,dispH);
    colormap(jet);
    if Dcomplex
        xlabel real;
        ylabel imaginary;
    else
        xlabel x;
        ylabel y;
    end
    
    
end

function [Ap]=SPLITsimple_simulate_pattern(S1, X, Y, T, alphas,f0, modp, Nph, showpattern)

Nang=length(alphas);
Z=Nph*Nang;
Ap=zeros(X,Y,Z);

%show pattern
for ang=1:Nang
    alpha=alphas(ang);
for k=1:Nph
    for j=1:X
        for i=1:Y
%     Ap(i,j,(ang-1)*Nph+k)=S1* (Pmin +  1*((cos(pi*( (j*cos(alpha)-i*sin(alpha)) /T + k/5 +f0(ang)  )))^2 ) )/(1+Pmin)  ;
Ap(i,j,(ang-1)*Nph+k)=S1* ( 1 - modp*((sin(pi*( (j*cos(alpha)-i*sin(alpha)) /T + k/5 +f0(ang)  )))^2 )  )   ;
%     Ap(i,j,(ang-1)*Nph+k)=S1* ( 1 - modp*((sin(pi*( (j*cos(alpha)-i*sin(alpha)) /T + k/5 +f0(ang)  )))^2 ) - modp2*((sin(pi*( (j*cos(alpha)-i*sin(alpha))*(2)/T + k/5 +f0(ang)+f2  )))^2 )  )   ;
        end
    end
end
end

if showpattern=='Y'
figure
colormap(hot)
subplot(2,2,1)
imagesc(Ap(:,:,1))
% subplot(2,2,2)
% imagesc(Ap(:,:,6))
% subplot(2,2,3)
% imagesc(Ap(:,:,11))
end





end

function [ Ap]=SPLITsimple_simulate_pattern_h(S1, X, Y, T, h,alphas,f0, modp, Nph, showpattern)

% w=4.5;  % if px=70nm then w should be about 4.5px for 315nm
% T=4.5*2;
% modp is the fringe contrast = [0,1]

Nang=length(alphas);
Z=Nph*Nang;
Ap=zeros(X,Y,Z);

%show pattern
for ang=1:Nang
    alpha=alphas(ang);
for k=1:Nph
    for j=1:X
        for i=1:Y
%     Ap(i,j,(ang-1)*Nph+k)=S1* (Pmin +  1*((cos(pi*( (j*cos(alpha)-i*sin(alpha)) /T + k/5 +f0(ang)  )))^2 ) )/(1+Pmin)      ;
    Ap(i,j,(ang-1)*Nph+k)=S1* ( 1 - modp*((sin(h*pi*( (j*cos(alpha)-i*sin(alpha)) /T + k/5 +f0(ang)  )))^2 ) )   ; %test 26 nov
        end
    end
end
end

if showpattern=='Y'
figure
colormap(hot)
subplot(2,2,1)
imagesc(Ap(:,:,1))
subplot(2,2,2)
imagesc(Ap(:,:,6))
subplot(2,2,3)
imagesc(Ap(:,:,11))
end



end

function gf_corr=SPLITsimple_CorrectPhasor(gf_nocorr, g_ref, B,h)
gf_corr=gf_nocorr(:,:,1);
Nangles=size(gf_nocorr,3);
X=size(gf_nocorr,1);
Y=size(gf_nocorr,2);

gfanglemod=zeros('like',gf_nocorr);
gfanglephase=zeros('like',gf_nocorr);
for anglesim=1:Nangles
     for ix=1:X
        for j=1:Y
%             if B(ix,j)>0
            if 1>0    
%        gfanglephase(ix,j,anglesim)=abs( mod( angle(gf_nocorr(ix,j,anglesim)) -  angle(g_ref(ix,j,anglesim)) , 2*pi) -pi );
        
       
       if h==1
        gfanglephase(ix,j,anglesim)=( mod( angle(gf_nocorr(ix,j,anglesim)) -  angle(g_ref(ix,j,anglesim)) , 2*pi) ); % test 26 nov
       if gfanglephase(ix,j,anglesim)>pi
           gfanglephase(ix,j,anglesim)=(2*pi)-gfanglephase(ix,j,anglesim);
       end                                                                                                          % test 26 nov
       else
           gfanglephase(ix,j,anglesim)=( mod( angle(gf_nocorr(ix,j,anglesim)) -  angle(g_ref(ix,j,anglesim)) , 2*pi) ); % test 26 nov
       if gfanglephase(ix,j,anglesim)>pi
           gfanglephase(ix,j,anglesim)=(2*pi)-gfanglephase(ix,j,anglesim);
       end                                                                                                          % test 26 nov
           
           
%           if abs( mod( angle(gf_nocorr(ix,j,anglesim)),pi) - mod( angle(g_ref(ix,j,anglesim)) , pi) ) < pi/2 
%               gfanglephase(ix,j,anglesim)= abs( mod( angle(gf_nocorr(ix,j,anglesim)),pi) - mod( angle(g_ref(ix,j,anglesim)) , pi) ) ;
%           else
%               gfanglephase(ix,j,anglesim)= pi- abs( mod( angle(gf_nocorr(ix,j,anglesim)),pi) - mod( angle(g_ref(ix,j,anglesim)) , pi) ) ;
%           end
          
%         gfanglephase(ix,j,anglesim)=abs(abs( mod( angle(gf_nocorr(ix,j,anglesim)) -  angle(g_ref(ix,j,anglesim)) , 2*pi) -pi )-pi);   % da correggere
%          gfanglephase(ix,j,anglesim)=( mod( angle(gf_nocorr(ix,j,anglesim)) -  angle(g_ref(ix,j,anglesim)) , 1*pi) ); % test 26 nov
%          if gfanglephase(ix,j,anglesim)>pi/2
%            gfanglephase(ix,j,anglesim)=(pi)-gfanglephase(ix,j,anglesim);
%           end 
%             if gfanglephase(ix,j,anglesim)>pi/2
%                 if gfanglephase(ix,j,anglesim)>3/2*pi
%                     gfanglephase(ix,j,anglesim)=(2*pi)-gfanglephase(ix,j,anglesim);
%                 elseif gfanglephase(ix,j,anglesim)>pi
%                     gfanglephase(ix,j,anglesim)=gfanglephase(ix,j,anglesim)-pi;
%                 else
%                     gfanglephase(ix,j,anglesim)=pi-gfanglephase(ix,j,anglesim);
%                 end
%                 
%             end 
       end
       
       
       gfanglemod(ix,j,anglesim) = abs( gf_nocorr(ix,j,anglesim) ) ;  % modulation is the same in a rotation    
            end
        end
     end
end

% 
% gfmod=mean(gfanglemod,3);
% gfphase=mean(gfanglephase,3);
% for ix=1:X
%     for j=1:Y
% %         if B(ix,j)>0
%          if 1>0
%     gf_corr(ix,j)=gfmod(ix,j).*exp(-gfphase(ix,j)*1i);
%         end
%     end
% end

%averaging vectors
gfcompl=mean(gfanglemod.*exp(-gfanglephase*1i),3);
gf_corr=gfcompl;
% weighted averaging vectors
% mods=sum(gfanglemod,3);
% for ix=1:X
%     for j=1:Y
%     gfweights(ix,j,:)=gfanglemod(ix,j,:)/mods(ix,j);
%     end
% end
% gfcompl=mean(gfweights.*gfanglemod.*exp(-gfanglephase*1i),3);
% gf_corr=gfcompl;

end

function [f,flin]=SPLITSIM_fraction_2comp(phaseimg_nothr, min_ph,delta_ph, mod1_nothr, max_mod, over)
nc=2;
[X,Y]=size(phaseimg_nothr);
over=1; % oversampling outside function

g_in=max_mod*cos(min_ph) ;      % 
s_in=0*max_mod*sin(min_ph) ;  % use only g axis
Pin=[g_in s_in];
g_out= max_mod*cos(min_ph+delta_ph) ;  %
s_out= 0*max_mod*sin(min_ph+delta_ph) ;  % use only g axis
Pout=[g_out s_out] ;
V0=Pout-Pin;
gnc(:,:,1)=mod1_nothr.*cos(phaseimg_nothr);
gnc(:,:,2)=mod1_nothr.*sin(phaseimg_nothr);

flin=zeros(over*X,over*Y,2);
f=zeros(over*X,over*Y,2);
klog=4;

for i=1:X
    for j=1:Y
        xv=squeeze(gnc(i,j,:))';
        flin(i,j,2)=dot((xv-Pin),V0)./(norm(V0)^2);
        flin(i,j,1)=1-flin(i,j,2);
        f(i,j,2)=1/(1+exp(- klog* (flin(i,j,2)-0.5) )) ;
        f(i,j,1)=1-f(i,j,2);
    end
end

end

function bool=iseven(x)

if mod(x,2) == 0
bool=1;
else
bool=0;
end
end


