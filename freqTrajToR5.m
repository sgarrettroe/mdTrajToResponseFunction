function [mean_w,S,Sf,t1,w1] = freqTrajToR5(p,freq,dipole)
%mimic what I am doing in C

n_levels = length(freq);

dt = p.dt; %0.01; %ps
nt = p.nt; %32;
ntint = p.ntint; %4;
nsteps = p.nsteps;%1000;
nprotons = p.nprotons;%1019;
nskip = p.nskip;%15;

t2_t4_pairs = [0 0];
order = 1;

%calculate the mean freq so we can go to the rotating frame
mean_w = mean(freq{1}(:));
disp(['mean freq = ' num2str(mean_w)]);

for i = 1:n_levels
    freq{i} = freq{i} - mean_w;
end


S  = zeros(1,nt);
Sf = zeros(1,2*nt);

sdw = zeros(1,nt);

t2 = t2_t4_pairs(1);
t4 = t2_t4_pairs(2);
nt2 = t2/ntint;
nt4 = t4/ntint;

ntraject = 3*nt + nt2 + nt4;
dwint = zeros(1,ntraject);
muint = zeros(1,ntraject);


if nsteps<ntraject*ntint, error('ERROR: trajectory is too short for response function calculation!'),end
nsamples = floor( (nsteps-ntraject*ntint)/nskip );
disp(sprintf('The number of samples for this set of delays is %6.4i\n',nsamples));

for iproton = 1:nprotons
    for isample = 1:nsamples
        if mod(isample+nsamples*(iproton-1),500)==0
            disp(['trajectory ' num2str(isample+nsamples*(iproton-1))]);
        end
        
        %linear spectroscopy
        if order >= 1
            order_calc = 1;
            n_levels = (order_calc + 1)/2;
            
            i_level = 1;
            offset = (isample-1)*nskip;
            dw = freq{i_level}(offset + (1:ntraject*ntint),iproton);
            mu = dipole{i_level}(offset + (1:ntraject*ntint),iproton);
            
            for i = 1:ntraject
                ind = (1:ntint)+(i-1)*ntint;
                dwint(i) = trapz(dw(ind));
                muint(i) = trapz(mu(ind));
            end
            
            %here is a literal translation of what is going on in the C
            %code, though is it right?
            sdw(1) = 0;
            sdw(2:end) = cumsum(dwint(1:nt-1));
            
            t1_ind = 1:nt;
            S = S + muint(1)*(muint(t1_ind).*exp(1i*sdw));
            
        end
        
    end
end

S = S./(isample+nsamples*(iproton-1));
zeropad = 2*nt;
Sf = fftshift(real(sgrsfft(S,zeropad)));
t1 = (0:(nt-1))*ntint*dt;
w1 = fftFreqAxis(t1,...
    'time_units','ps',...
    'freq_units','cm-1',...
    'fftshift','on',...
    'zeropad',zeropad);
