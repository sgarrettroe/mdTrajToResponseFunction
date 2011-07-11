function d = loadR5(parameter_file_name,time_file_name,fname)

i_pix = 16;
range = [-400 400];
flag_debug = 1;

p = read_parameter_file(parameter_file_name);

%initialize 3d data structures
d = load3d;

nt = p.nt;
ntint = p.ntint;
d.nt = p.nt;
dt = p.dt*p.ntint;
n_steps = p.nsteps_in_file; %this is going to break!
d.time = (0:nt-1)*dt;
d.freq = []; %this means no spectrometer
d.freq_units = 'cm-1';
d.time_units = 'ps';

dummy = load(time_file_name);
t2 = dummy(1);
t4 = dummy(2);

d.t2 = t2; %ps!
d.t4 = t4;

d.phase = 0;
d.undersampling = 0;
d.basename = fname;
d.comment = [...
    ' nt = ' num2str(nt) ...
    ' ntint = ' num2str(ntint) ...
    ' dt = ' num2str(dt) ...
    ' nsteps = ' num2str(n_steps)];


%initialize the 
R1 = zeros(nt,nt,nt);
R2 = zeros(nt,nt,nt);
R3 = zeros(nt,nt,nt);
R4 = zeros(nt,nt,nt);

%default zeropadded length
n_zp = 2*nt;
fft_type = 'sgrsifft';
apodization = 'none';

dummy = zeros(nt^3,8); %8 columns for 4 diagrams real and imag
if exist(fname,'file')
    try
        dummy = dummy + load(fname);
    catch
        warning(['WARNING!!! Something is wrong with ', fname]);
    end
else
    warning(['WARNING!!! A file is missing! I can''t find ', fname]);
end

count =0;
for j=1:nt;
    for k=1:nt;
        for l=1:nt;
            count=count+1;
            R1(k,j,l) = dummy(count,1) + 1i*dummy(count,2);
            R2(k,j,l) = dummy(count,3) + 1i*dummy(count,4);
            R3(k,j,l) = dummy(count,5) + 1i*dummy(count,6);
            R4(k,j,l) = dummy(count,7) + 1i*dummy(count,8);
        end
    end
end
d.R1 = R1;
d.R2 = R2;
d.R3 = R3;
d.R4 = R4;

%fix frequency axes for shift_w
switch p.n_levels
    case 1
        shift_w3 = 0;
        shift_w5 = 0;
    case 2
        shift_w = (p.mean_w(1)-p.mean_w(2))/2;
        shift_w3 = -shift_w;
        shift_w5 = 0;
    case 3
        shift_w = (p.mean_w(1)-p.mean_w(2))/2;
        mw = mean(p.mean_w);
        shift_w3 = -shift_w;
        shift_w5 = -(p.mean_w(1) - mw);
end

d = absorptive3d(d,d.phase,i_pix,n_zp,fft_type,apodization,range,'on');

d.w1 = d.w1 + p.mean_w(1);
d.w3 = d.w3 + p.mean_w(1) + shift_w3;
d.w5 = d.w5 + p.mean_w(1) + shift_w5;

if flag_debug,figure(1),my3dPlot(d.w1,d.w3,d.w5,d.R);end

  %the spectrum is totally inverted, so fix that
  %the lines above don't do the trick so we do it now
  %I don't understand why this is necessary...
  %d(i).R = flipdim(flipdim(flipdim(d(i).R,1),2),3);
  %d.w1 = fliplr(d.w1);
  %d.w3 = fliplr(d.w3);
  %d.w5 = fliplr(d.w5);
  %if flag_debug,figure(2),my3dPlot(d.w1,d.w3,d.w5,d.R);end
end
