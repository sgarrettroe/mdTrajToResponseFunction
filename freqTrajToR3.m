function [s1d,s2d] = freqTrajToR3(p,f,x)
global c 

n_contours = 10;
flag_noncondon = p.flag_noncondon;
flag_bsxfun = exist('bsxfun','builtin'),
flag_plot = 1;
  
tic
% try again
%nt = 20;
%nskip = 1;
%ntint = 4;
nt = p.nt;
nskip = p.nskip;
ntint = p.ntint;
nsteps = p.nsteps;
order = p.order;
t2 = p.t2;
flag_twolevelsystem = p.flag_twolevelsystem;


freq = f{1};
if isempty(x)
    mean_mu = 1;
else
    if iscell(x)
        mu = x{1};
        mean_mu = mean(mu);
    else
        if length(x)==1
            flag_noncondon = 0;
            mean_mu = x;
        elseif length(x)==length(freq)
            mu = x;
            mean_mu = mean(mu);
        end
    end
end

%how many independent microtrajectories from this total trajectory
ntrajectories = floor((nsteps/ntint - nt)/nskip);

%do the calculation in the rotating frame of the 01 transition
mean_w = mean(f{1});
dt = 0.01*ntint;

%convert from wavenumbers to rad/ps-1
dw = (freq-mean_w)*2*pi*c*1e-12;

%take every 4th point
%dw = dw(1:ntint:end);

%take an average over ntint points
dw = dw(1:floor(nsteps/ntint)*ntint); %throw away extra points
nsteps = length(dw); %update nsteps
dw = trapz(reshape(dw,ntint,[]),1)./(ntint-1); 

mu = mu(1:floor(nsteps/ntint)*ntint); %throw away extra points
mu = trapz(reshape(mu,ntint,[]),1)./(ntint-1); 

%1D spectrum
S  = zeros(1,nt);
Sf = zeros(1,2*nt);

%build a matrix DW1 from the frequency trajectory consisting of each block of
%that trajectory that we will integrate and average. Each column has nt
%points with a different time origin given by <it0>
DW1 = zeros(nt,ntrajectories);
it0 = 1:nskip:nsteps/ntint-nt;
for is = 1:ntrajectories
    DW1(:,is) = dw(it0(is):it0(is)+nt-1);
end

%calculate the time integral for each minitrajectory
sdw = cumtrapz(DW1,1)*dt;

%then calcuate the ensemble average of exp(1i <time integral>)
if flag_noncondon==0
    disp('condon linear');
    S = mean_mu^2.*sum(exp(-1i*sdw),2);
else
    disp('non-condon linear');

    %try bsxfun instead???
    if flag_bsxfun == 0
      MU0 = repmat(mu(it0),[nt,1]);
    end

    MU1 = zeros(nt,ntrajectories);
    for is = 1:ntrajectories
        MU1(:,is) = mu(it0(is):it0(is)+nt-1);
    end
    if flag_bsxfun == 0
      %keyboard
      S = sum(MU1.*MU0.*exp(-1i*sdw),2);
    else
      S = sum(bsxfun(@times,mu(it0),MU1).*exp(-1i*sdw),2);
    end
end

%normalize
S = S*nskip/(nsteps/ntint-nt);

%the first time points need to be divided by 2 (Sect 9.5.3)
S(1) = S(1)/2;

%do the fft
Sf = real(fftshift(fft(S,2*nt)));

%plot
t1 = (0:nt-1)*dt;
t3 = t1;
w1 = fftFreqAxis(t1,...
    'time_units','ps',...
    'freq_units','cm-1',...
    'zeropad',2*nt,...
    'fftshift','on');
w1 = w1 + mean_w;
s1d = construct1d;
s1d.w1 = w1 ;
s1d.t1 = t1;
s1d.R1 = S;
s1d.R  = Sf;
s1d.centerfreq = mean_w;


if (order<3)
    s2d = construct2d;
    return;
end

clear dw DW1 MU1 mu sdw

%2D spectrum

%do the calculation in the rotating frame for each transition and
%then add the anharmonic frequency shifts at the end
mean_w01 = mean(f{1});
mean_w12 = mean(f{2});
Delta_anh = mean_w01 - mean_w12;

if flag_twolevelsystem
  %if only the 10 peak is calculated stay centered on it
  mean_w = mean_w01;
  shift_w = 0;
  Delta_anh = 0;
else
  %but if we are doing multilevel then stay centered between 01
  %and 12 transitions
  mean_w = (mean_w01 + mean_w12)/2;
  shift_w = Delta_anh/2;
end


%convert from wavenumbers to rad/ps-1
dw_01 = (f{1}-mean_w01)*2*pi*c*1e-12;
dw_12 = (f{2}-mean_w12)*2*pi*c*1e-12;
shift_w = shift_w*2*pi*c*1e-12;

%set up the time axes
[T1,T3] = meshgrid(t1,t1);

%take an average over ntint points
dw_01 = dw_01(1:floor(nsteps/ntint)*ntint); %throw away extra points
dw_12 = dw_12(1:floor(nsteps/ntint)*ntint); %throw away extra points
dw_01 = trapz(reshape(dw_01,ntint,[]),1)./(ntint-1); 
dw_12 = trapz(reshape(dw_12,ntint,[]),1)./(ntint-1); 

%take an average over ntint points
mu_01 = x{1}(1:floor(nsteps/ntint)*ntint); %throw away extra points
mu_12 = x{2}(1:floor(nsteps/ntint)*ntint); %throw away extra points
mu_01 = trapz(reshape(mu_01,ntint,[]),1)./(ntint-1); 
mu_12 = trapz(reshape(mu_12,ntint,[]),1)./(ntint-1); 

nsteps = length(dw_01); %update nsteps


%
it2 = round(t2/dt);

ntrajectories = floor((nsteps/ntint - 2*nt - it2)/nskip);
R_r  = zeros(nt,nt); % rephasing diagrams (really R1 + R2)
R_nr = zeros(nt,nt); % non-rephasing (R4 + R5)
R3  = zeros(nt,nt); % rephasing excited state abs (R3)
R6  = zeros(nt,nt); % non-rephasing e.s.a. (R6)
R = zeros(2*nt,2*nt); % final fourier transformed spectrum

%build a matrix DW1 from the frequency trajectory consisting of each block of
%that trajectory that we will integrate and average. Each column has nt
%points with a different time origin given by <it0>. For DW2 we do the
%same but as a function of t1 and t3
DW1 = zeros(nt,ntrajectories);
DW2 = zeros(nt,nt,ntrajectories);

it0 = 1:nskip:nsteps/ntint-nt;
for is = 1:ntrajectories
    DW1(:,is) = dw_01(it0(is):it0(is)+nt-1);
end
for is = 1:ntrajectories
    for it1 = 1:nt
        DW2(it1,:,is) = dw_01(it0(is)+it1+it2:it0(is)+it1+it2+nt-1);
    end
end

%calculate the time integral for each minitrajectory
sdw = cumtrapz(DW1,1)*dt;

%copy this out into a matrix for (t1,t3)
if flag_bsxfun == 0
  SDW = repmat(reshape(sdw,[nt 1 ntrajectories]),[1 nt 1]);
else
  SDW = reshape(sdw,[nt 1 ntrajectories]);
end

%for the t3 dependence
%PDW = zeros(nt,nt,ntrajectories);
PDW = cumtrapz(DW2,2)*dt;

if flag_noncondon == 0
    disp('condon ground state + stim emission');
    %The factor of 2 is for the 2 population pathways. The bsxfun applies
    %the operation (plus, minus or times) and automatically expands the
    %singleton dimensions of the smaller matrix to save memory (ie you
    %don't need to do a repmat to make a SDW an nt x nt matrix and the same
    %for MU1_MU0)    
    if flag_bsxfun == 0
      R_r  = 2*mean_mu^4.*sum(exp(-1i.*PDW + 1i*SDW),3);
      R_nr = 2*mean_mu^4.*sum(exp(-1i.*PDW - 1i*SDW),3);
    else
      R_r  = 2*mean_mu^4.*sum(exp(-1i.*(bsxfun(@minus,PDW, SDW))),3);
      R_nr = 2*mean_mu^4.*sum(exp(-1i.*(bsxfun(@plus, PDW, SDW))),3);
    end
else
    disp('non-condon ground state + stim emission');
    
    MU1_MU0 = zeros(nt,ntrajectories);
    for is = 1:ntrajectories
        MU1_MU0(:,is) = mu_01(it0(is)).*mu_01(it0(is):it0(is)+nt-1);
    end
    if flag_bsxfun == 0
      MU1_MU0 = repmat(reshape(MU1_MU0,[nt 1 ntrajectories]),[1 nt 1]);
    else
      MU1_MU0 = reshape(MU1_MU0,[nt 1 ntrajectories]);
    end
    
    MU3_MU2 = zeros(nt,nt,ntrajectories);
    for is = 1:ntrajectories
        for it1 = 1:nt
            MU3_MU2(it1,:,is) = mu_01(it0(is)+it1+it2).*mu_01(it0(is)+it1+it2:it0(is)+it1+it2+nt-1);
        end
    end    
    
    %The factor of 2 is for the 2 population pathways. The bsxfun applies
    %the operation (plus, minus or times) and automatically expands the
    %singleton dimensions of the smaller matrix to save memory (ie you
    %don't need to do a repmat to make a SDW an nt x nt matrix and the same
    %for MU1_MU0)
    if flag_bsxfun == 0
      R_r  = 2*sum(MU3_MU2.*MU1_MU0.*exp(-1i.*PDW + 1i.*SDW),3);
      R_nr = 2*sum(MU3_MU2.*MU1_MU0.*exp(-1i.*PDW - 1i.*SDW),3);
    else
      R_r  = 2*sum(bsxfun(@times,MU3_MU2,MU1_MU0).*exp(-1i.*(bsxfun(@minus,PDW, SDW))),3);
      R_nr = 2*sum(bsxfun(@times,MU3_MU2,MU1_MU0).*exp(-1i.*(bsxfun(@plus, PDW, SDW))),3);
    end
end    
R_r  = exp( 1i*shift_w*T3).*R_r;
R_nr = exp( 1i*shift_w*T3).*R_nr;

if flag_twolevelsystem == false
    %add excited state if we need to

    %DW1 is the same from the 01 calc above
    
    %DW2 is now the 12 transition
    for is = 1:ntrajectories
        for it1 = 1:nt
            DW2(it1,:,is) = dw_12(it0(is)+it1+it2:it0(is)+it1+it2+nt-1);
        end
    end
    
    %for the t3 dependence
    %PDW = zeros(nt,nt,ntrajectories);
    PDW = cumtrapz(DW2,2)*dt;

   if flag_noncondon == 0
       disp('condon excited state abs');
       if flag_bsxfun == 0
	 R3 = -2 * mean_mu^4 .* sum(exp(-1i.*PDW + 1i*SDW),3);
	 R6 = -2 * mean_mu^4 .* sum(exp(-1i.*PDW - 1i*SDW),3);
       else
	 R3 = -2 * mean_mu^4 .* sum(exp(-1i.*(bsxfun(@minus,PDW, SDW))),3);
	 R6 = -2 * mean_mu^4 .* sum(exp(-1i.*(bsxfun(@plus, PDW, SDW))),3);
       end
       %R_r  = sum(exp(-1i.*PDW + 1i*SDW),3);
       %R_nr = sum(exp(-1i.*PDW - 1i*SDW),3);
       %R_r  = R_r  - mean_mu^4.*sum(exp(-1i.*(bsxfun(@minus,PDW, SDW))),3);
       %R_nr = R_nr - mean_mu^4.*sum(exp(-1i.*(bsxfun(@plus, PDW, SDW))),3);
   else
       disp('non-condon excited state abs');
       
       MU1_MU0 = zeros(nt,ntrajectories);
       for is = 1:ntrajectories
           MU1_MU0(:,is) = mu_12(it0(is)) .* mu_12(it0(is):it0(is)+nt-1);
       end
       if flag_bsxfun == 0
	 MU1_MU0 = repmat(reshape(MU1_MU0,[nt 1 ntrajectories]),[1 nt 1]);
       else
	 MU1_MU0 = reshape(MU1_MU0,[nt 1 ntrajectories]);
       end

       
       MU3_MU2 = zeros(nt,nt,ntrajectories);
       for is = 1:ntrajectories
           for it1 = 1:nt
               MU3_MU2(it1,:,is) = mu_12(it0(is)+it1+it2).*mu_12(it0(is)+it1+it2:it0(is)+it1+it2+nt-1);
           end
       end
    
       %The factor of 2 is not needed!  it is handled by the mu. The bsxfun applies
       %the operation (plus, minus or times) and automatically expands the
       %singleton dimensions of the smaller matrix to save memory (ie you
       %don't need to do a repmat to make a SDW an nt x nt matrix and the same
       %for MU1_MU0)
       if flag_bsxfun == 0
	 R3 = -1 * sum(MU3_MU2.*MU1_MU0.*exp(-1i.*PDW + 1i.*SDW),3);
	 R6 = -1 * sum(MU3_MU2.*MU1_MU0.*exp(-1i.*PDW - 1i.*SDW),3);
       else
	 R3 = -1 * sum(bsxfun(@times,MU3_MU2,MU1_MU0).*exp(-1i.*(bsxfun(@minus,PDW, SDW))),3);
	 R6 = -1 * sum(bsxfun(@times,MU3_MU2,MU1_MU0).*exp(-1i.*(bsxfun(@plus, PDW, SDW))),3);
       end
   end
   R3 = exp(-1i*shift_w.*T3).*R3;
   R6 = exp(-1i*shift_w.*T3).*R6;

end
  
%normalize
R_r  = R_r *nskip/(nsteps/ntint-nt);
R_nr = R_nr*nskip/(nsteps/ntint-nt);
R3   = R3  *nskip/(nsteps/ntint-nt);
R6   = R6  *nskip/(nsteps/ntint-nt);

%the first time points need to be divided by 2 (Sect 9.5.3)
%R_r(:,1)  = R_r(:,1)/2;
%R_r(1,:)  = R_r(1,:)/2;
%R_nr(:,1) = R_nr(:,1)/2;
%R_nr(1,:) = R_nr(1,:)/2;

%do the fft
R_nr_f = fftshift(sgrsfft2(R_nr + R6,2*nt));
R_r_f  = fftshift(sgrsfft2(R_r  + R3,2*nt));
R_r_f  = fliplr(  circshift(R_r_f,[0 -1]));
%R_r_f  = circshift(fliplr(R_r_f),[-1 0]);
R      = real(R_r_f + R_nr_f);

%the freq axis for plotting
w1 = fftFreqAxis(t1,...
    'time_units','ps',...
    'freq_units','cm-1',...
    'zeropad',2*nt,...
    'fftshift','on');
w1 = w1 + mean_w01 ;
w3 = w1 - Delta_anh/2;

s2d = construct2d;
s2d.t1 = t1;
s2d.t3 = t1;
s2d.w1 = w1;
s2d.w3 = w3;
s2d.t2 = it2;
s2d.R1 = R_nr;
s2d.R2 = R_r;
s2d.R3 = R3;
s2d.R6 = R6;
s2d.R  = R;
s2d.zeropad = 2*nt;


%
toc
disp('done')

if flag_plot
%figure(10)
%contourf(real(R_r))
%figure(11)
%contourf(real(R_nr))

%
%add the excited state absorption
% [T1,T3] = meshgrid(t1,t3);
% Delta = 240;
% R_r_anh = (2*exp(1i*Delta/2*(2*pi*c*1e-12).*T3) - 2*exp(-1i*Delta/2/(2*pi*c*1e-12).*T3)).*R_r;
% R_nr_anh = (2*exp(1i*Delta/2*(2*pi*c*1e-12).*T3) - 2*exp(-1i*Delta/2/(2*pi*c*1e-12).*T3)).*R_nr;

% %the first time points need to be divided by 2 (Sect 9.5.3)
% R_r_anh(:,1) = R_r_anh(:,1)/2;
% R_r_anh(1,:) = R_r_anh(1,:)/2;
% R_nr_anh(:,1) = R_nr_anh(:,1)/2;
% R_nr_anh(1,:) = R_nr_anh(1,:)/2;
% 
% %do the fft
% R_nr_anh_f = fftshift(fft2(R_nr_anh,2*nt,2*nt));
% R_r_anh_f  = fftshift(fft2(R_r_anh, 2*nt,2*nt));
% R_r_anh_f = fliplr(circshift(R_r_anh_f,[0 -1]));
% %R_r_anh_f  = circshift(fliplr(R_r_anh_f),[-1 0]);
% R_anh    = real(R_r_anh_f + R_nr_anh_f);

%figure(12)
%contourf(w1,w1,R,n_contours)

%figure(1),
%plot((0:nt-1)*ntint,dw(1:nt),'o',...
%    (0:nt-1)*ntint,DW1(:,1),'s',...
%    0:(nt-1)*ntint,(freq(1:(nt-1)*ntint+1)-mean_w)*2*pi*c*1e-12);

figure(10)
title('FID');
plot(t1,real(S),'-o')
xlabel('t1')
ylabel('S(t)')

figure(11)
title('1D spectrum');
plot(w1,Sf,'-o')
xlabel('\omega_1')
ylabel('S(\omega)')

%
figure(12),clf
my2dPlot(w1,w3,R,'n_contours',n_contours,'pumpprobe',false)

if flag_twolevelsystem==0
  figure(13),clf
  R_nr_f = fftshift(sgrsfft2(R_nr,2*nt));
  R_r_f  = fftshift(sgrsfft2(R_r,2*nt));
  R_r_f  = fliplr(  circshift(R_r_f,[0 -1]) );
  R = real(R_r_f + R_nr_f);
  my2dPlot(w1,w3,R,'n_contours',n_contours,'pumpprobe',false)

  figure(14),clf
  R_nr_f = fftshift(sgrsfft2(R6,2*nt));
  R_r_f  = fftshift(sgrsfft2(R3,2*nt));
  R_r_f  = fliplr(  circshift(R_r_f,[0 -1]) );
  R = real(R_r_f + R_nr_f);
  my2dPlot(w1,w3,R,'n_contours',n_contours,'pumpprobe',false)
end

end %if flag_plot