%-- Function: calc_a_b.m
%-- Inputs:  freq (vector NF x 1)
%--             dir (ND x 1)
%--             energy2d:  (time, freq, dir)
%-------------------------------------------------------------------------
function [ data ] = WIS_calc_ab(wisdirs,freqs,energy2d)

data = {};

rdir = deg2rad(wisdirs);
dtheta = abs(rdir(2)-rdir(1));


%-- input energy is 3 Dims (time, freq, dir)
ds = energy2d;   % (ND x NF)
[ND, NF, NT ] = size(ds);   %-- # Time, # Freqs, # Dir
bands = NF;
dbins = ND;

deg_per_bin = 360/dbins;

bw = diff(freqs);
bw(end+1) = bw(end);

% %% Assign WW3 directions to correct bins (direction coming from, starting with 5 deg bin)
% angle = round(rad2deg(rdir)) + 180;
% idx = find(angle >= 360);
% angle(idx) = angle(idx) - 360;
% dir_idx = round(angle./deg_per_bin);
% bin_dirs = zeros(dbins,1);
% for i = 1:length(dir_idx)
%     idx = dir_idx(i);
%     bin_dirs(idx) = angle(i);
% end
bin_dirs = wisdirs;


%% Initialize a,b params
icnt=0;
pstnid = '';
a0 = zeros(NF,NT);
a1 = zeros(NF,NT);
b1 = zeros(NF,NT);
a2 = zeros(NF,NT);
b2 = zeros(NF,NT);
energy = zeros(NF,NT);
Hs = [];
Tp = [];
Dp = [];
Ta = [];

        
zero_energy = [];
energy = zeros(NF,1);

sp1d = zeros(NF,NT);
ds_set = zeros(dbins,bands,NT);

for icnt = 1:NT
    ds = squeeze(energy2d(:,:,icnt));       %-- (ND x NF)
              
    for i = 1:bands
        for j = 1:dbins
            ds_set(j,i,icnt) = ds(j,i);
        end
    end   
    
    
    %--   When calculating moments, rotate directions by pi so that the resulting 
    %-   coefficients are in true compass "arriving from" coordinates.
    zero_energy(icnt) = true;
    for i = 1:bands
        
        for j = 1:dbins
            
            %--   Little floating point exception check here to catch 
            %--   very small double-precision energies from ww3 model
            if(ds_set(j,i,icnt) < 1.0e-15) 
                ds_set(j,i,icnt) = 0.0;
            end
            
            
            a0(i,icnt) = a0(i,icnt)+ds_set(j,i,icnt);
            a1(i,icnt) = a1(i,icnt)+ds_set(j,i,icnt)*cos(deg2rad(bin_dirs(j)));
            b1(i,icnt) = b1(i,icnt)+ds_set(j,i,icnt)*sin(deg2rad(bin_dirs(j)));
            a2(i,icnt) = a2(i,icnt)+ds_set(j,i,icnt)*cos(deg2rad(2.*bin_dirs(j)));
            b2(i,icnt) = b2(i,icnt)+ds_set(j,i,icnt)*sin(deg2rad(2.*bin_dirs(j)));
            ds_set(j,i,icnt) = ds_set(j,i,icnt) / deg_per_bin;            
                 
        end     %-- End of 1st loop through directions (dbins)
        
        %--   Normalize fourier coeffiecients
        if(a0(i,icnt) > 0) 
            zero_energy(icnt) = false;
            a1(i,icnt) = a1(i,icnt)/a0(i,icnt);
            b1(i,icnt) = b1(i,icnt)/a0(i,icnt);
            a2(i,icnt) = a2(i,icnt)/a0(i,icnt);
            b2(i,icnt) = b2(i,icnt)/a0(i,icnt);
            wisdirs(i,icnt) = rad2deg(atan2(b1(i,icnt),a1(i,icnt)));
            if (wisdirs(i,icnt) < 0) 
                wisdirs(i,icnt) = wisdirs(i,icnt) + 360;
            end
        end     
        data.a0(i,icnt) = a0(i,icnt);
        data.a1(i,icnt) = a1(i,icnt);
        data.b1(i,icnt) = b1(i,icnt);
        data.a2(i,icnt) = a2(i,icnt);
        data.b2(i,icnt) = b2(i,icnt);
        data.energy(i,icnt) = sum(ds(:,i))*dtheta;
        
        
    end  %- End of freq loop
    
    %--   Calculate Hs, Tp, Dp, Ta

    if (~ zero_energy(icnt)) 
        total_energy = 0;
        M0 = 0;
        M1 = 0;
        for i = 1:bands
            energy(i) = a0(i,icnt) * bw(i);
            total_energy = total_energy + energy(i);
            M0 = M0 + energy(i);
            M1 = M1 + energy(i) * freqs(i);
        end
        Hs(icnt) = 4 * sqrt(total_energy);
        [tmp,peak_band]=max(a0(:,icnt));
        Tp(icnt) = 1.0 / freqs(peak_band);
        Dp(icnt) = wisdirs(peak_band,icnt);
        Ta(icnt) = M0 / M1;
       
    end
    
end   
data.hs = Hs;
data.tp = Tp;
data.dp = Dp;
data.ta = Ta; 
end