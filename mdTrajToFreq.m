function [freq,mu,field] = mdTrajToFreq(p,force_file_name,coord_file_name,varargin)
%mimic what I am doing in C

n_levels = 2;
flag_massweightedforces = p.flag_massweightedforces;

%nsteps = 100;
%nmols = 1019;
%q_H = 0.4238;

nsteps = p.nsteps;
nmols = p.nmols;
q_H = p.q_H;

%Skinner & Corcelli give units of field in au, which are 
%(field) 1 au = E_h e^-1 a_0^-1 (Hartree/(charge.bohr))
%
% but gromacs is in kJ/mol so the field is
%(field) 1 gmx = 1 kJ mol-1 e-1 nm-1
%
% so I get a conversion factor of 
%
% E_h      4.359744e-18 J    1 kJ     6.022e23
% ----  =  --------------- -------- ------------- = 49613.323 kJ / mol e nm
% e a_0     e 0.052918 nm   1000 J      1 mol

a = p.a;
b = p.b;
mu_mug = p.mu_mug;

force_fid = fopen(force_file_name);
coord_fid = fopen(coord_file_name);

posO = zeros(3,1);
posH1 = zeros(3,1);
posH2 = zeros(3,1);
forceO = zeros(3,1);
forceH1 = zeros(3,1);
forceH2 = zeros(3,1);

freq = cell(1,n_levels);
x    = cell(1,n_levels);
mu   = cell(1,n_levels);

field = zeros(nsteps,nmols);
for i = 1:n_levels
    freq{i} = zeros(nsteps,nmols);
    mu{i}   = zeros(nsteps,nmols);
end

%for mass weighted forces
mH = 2;
mO = 16;
red_mass = mH*mO/(mH+mO);

for i = 1:nsteps
    if mod(i,1000)==0,disp(['step ' num2str(i)]);end
    t = fscanf(coord_fid,'%f',1);
    t = fscanf(force_fid,'%f',1);
    for j = 1:nmols
        posO  = fscanf(coord_fid,'%f',3);
        posH1 = fscanf(coord_fid,'%f',3);
        posH2 = fscanf(coord_fid,'%f',3);

        forceO  = fscanf(force_fid,'%f',3);
        forceH1 = fscanf(force_fid,'%f',3);
        forceH2 = fscanf(force_fid,'%f',3);
        
        bond = posH1-posO;
        bond_length = sqrt(bond'*bond);
        u_bond = bond./bond_length;
        
        if flag_massweightedforces
            force = (forceH1/mH - forceO/mO)*red_mass;
        else
            %forces on H1
            force = forceH1;
        end
        
        %projected on the OH bond vector
        force_proj = force'*u_bond;
        
        %convert force to field
        field(i,j) = force_proj/q_H;
        
        %convert field to frequency w_01
        freq{1}(i,j) = a{1}(1) + a{1}(2).*field(i,j);

        %convert field to frequency w_12
        freq{2}(i,j) = a{2}(1) + a{2}(2).*field(i,j);
        
        %non-Condon effects
        x{1} = (b{1}(1) + b{1}(2)*freq{1}(i,j));
        x{2} = (b{2}(1) + b{2}(2)*freq{2}(i,j));
        %x{1} = (freq{1}(i,j)<2500);
        %x{2} = (freq{1}(i,j)<2500);
        mu{1}(i,j) = (mu_mug(1) + mu_mug(2)*field(i,j))*x{1};
        mu{2}(i,j) = (mu_mug(1) + mu_mug(2)*field(i,j))*x{2};
        
        
    end
end

