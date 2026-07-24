function [rad]=fRadTranCanopy_Matrix_tridiagonal_with_wood_for_release(sunang,rad,waveband,prm)
% Supplemental program 14.3 from Bonan 2019 'Climate change and
% terrestrial ecosystem modeling'

% --------------------------------------------
% Calculate and graph light profiles in canopy
% --------------------------------------------
%clear
% --- Model parameters

params.numrad = 1;         % Number of wavebands (visible, near-infrared)
params.vis = 1;            % Array index for visible waveband
params.nir = 2;            % Array index for near-infrared waveband
params.sun = 1;            % Array index for sunlit leaf
params.sha = 2;            % Array index for shaded leaf
params.npts = size(sunang.theta_rad,1);           % Number of grid points to process

%prm.dfw is wood area index profile
%prm.markov_wood is omega wood profile


switch waveband
    case 'par'

        rho = prm.par_reflect;     %: MB this params.vis does not work anymore               % Leaf reflectance (visible)

        tau = prm.par_trans;                    % Leaf transmittance (visible)

        rho_wood=prm.par_reflect_wood;

        omega = 1-prm.par_absorbed;      % Leaf scattering coefficient for canopy , MB: careful, omega here is not clumping
        flux.albsoib = prm.par_soil_refl;                         % Direct beam albedo of ground (visible)
        flux.albsoid = flux.albsoib;

    case 'nir'
        rho = prm.nir_reflect;     %: MB this params.vis does not work anymore               % Leaf reflectance (visible)
        tau = prm.nir_trans;                    % Leaf transmittance (visible)
        omega = 1-prm.nir_absorbed;
        rho_wood=prm.nir_reflect_wood;
        flux.albsoib = prm.nir_soil_refl;                         % Direct beam albedo of ground (visible)
        flux.albsoid = flux.albsoib;
end


% --- Model options

light = 'Norman';        % Use Norman radiative transfer
% light = 'Goudriaan';     % Use Goudriaan radiative transfer
% light = 'TwoStream';     % Use two-stream approximation radiative transfer


canopy.nsoi = 1;                                     % First layer is soil
canopy.nbot = canopy.nsoi + 1;                    % Bottom leaf layer

% --- Define plant canopy


% Set canopy LAI, layer LAI increment, and number of layers


LAI_profile=prm.dff;
WAI_profile=prm.dfw;
canopy.nveg=length(LAI_profile);
canopy.ntop = canopy.nbot + canopy.nveg - 1;   % Top leaf layer
for i=canopy.ntop:-1:canopy.nbot
    canopy.dlai(i)=LAI_profile(i-1); %MB: pushing vector up by one so that soil is layer 1
    canopy.dwai(i)=WAI_profile(i-1);
    canopy.clumpfac(i)=prm.markov(i-1);
    canopy.clumpfac_wood(i)=prm.markov_wood(i-1);
    leafanglecode(i)=prm.leafanglecode(i-1);
    leafanglecode_wood(i)=prm.leafanglecode_wood(i-1);

end

canopy.lai=sum(LAI_profile);


for i=canopy.nbot:length(canopy.dlai) %***********MB: LAI profile has now layer 1 as soil, so LAI layer 1 is 0
    switch leafanglecode(i)
        case 1
            %prm.leafanglecode(i)=1;
            %g(1,:)=[0.217993772,  0.205566574, 0.181747777, 0.149126381, 0.111758486, 0.07430429,  0.040992495, 0.016051798, 0.002458428];
            g(1,:)=[0.220605996,0.207341712,0.182558443,0.149048672,0.111175648,0.073128095,0.039618324,0.014835055,0.001688055];
        case 2
            %prm.leafanglecode(i)=2;
            %g(2,:)=[0.013816774, 0.044002464, 0.072956086, 0.099621579, 0.1232069, 0.143184018, 0.158760891, 0.169409487, 0.175041802];
            g(2,:)=[0.01516472,0.045145546,0.073731916,0.100052292,0.12323514,0.142757539,0.157922259,0.168380687,0.173609901];
        case 3
            %prm.leafanglecode(i)=3;
            %g(3,:)=[0.001813, 0.014469, 0.038908, 0.072257, 0.11037, 0.148837, 0.182804, 0.208301, 0.222241];
            g(3,:)=[0.001688055,0.014835055,0.039618324,0.073128095,0.111175648,0.149048672,0.182558443,0.207341712,0.220605996];
        case 4
            %prm.leafanglecode(i)=4;
            g(4,:)=[0.007084426, 0.053135816, 0.126252094, 0.19247555,  0.220395847, 0.197361602, 0.133842925, 0.059766886, 0.009684853];
        case 5
            %prm.leafanglecode(i)=5;

        case 6
            %prm.leafanglecode(i)=6;
            g(6,:)=[0.111111111,  0.111111111, 0.111111111, 0.111111111, 0.111111111, 0.111111111, 0.111111111, 0.111111111, 0.111111111];
    end
    switch leafanglecode_wood(i)
        case 1
            %prm.leafanglecode(i)=1;
            %g(1,:)=[0.217993772,  0.205566574, 0.181747777, 0.149126381, 0.111758486, 0.07430429,  0.040992495, 0.016051798, 0.002458428];
            g_wood(1,:)=[0.220605996,0.207341712,0.182558443,0.149048672,0.111175648,0.073128095,0.039618324,0.014835055,0.001688055];
        case 2
            %prm.leafanglecode(i)=2;
            %g(2,:)=[0.013816774, 0.044002464, 0.072956086, 0.099621579, 0.1232069, 0.143184018, 0.158760891, 0.169409487, 0.175041802];
            g_wood(2,:)=[0.01516472,0.045145546,0.073731916,0.100052292,0.12323514,0.142757539,0.157922259,0.168380687,0.173609901];
        case 3
            %prm.leafanglecode(i)=3;
            %g(3,:)=[0.001813, 0.014469, 0.038908, 0.072257, 0.11037, 0.148837, 0.182804, 0.208301, 0.222241];
            g_wood(3,:)=[0.001688055,0.014835055,0.039618324,0.073128095,0.111175648,0.149048672,0.182558443,0.207341712,0.220605996];
        case 4
            %prm.leafanglecode(i)=4;
            g_wood(4,:)=[0.007084426, 0.053135816, 0.126252094, 0.19247555,  0.220395847, 0.197361602, 0.133842925, 0.059766886, 0.009684853];
        case 5
            %prm.leafanglecode(i)=5;

        case 6
            %prm.leafanglecode(i)=6;
            g_wood(6,:)=[0.111111111,  0.111111111, 0.111111111, 0.111111111, 0.111111111, 0.111111111, 0.111111111, 0.111111111, 0.111111111];
    end
end




theta_l=deg2rad([5, 15, 25, 35, 45, 55, 65, 75, 85]);



%Feb 2025, following page 102 research book#5 after Smolander et al 2003,
%modify leaf optical properties to accoutn for increased recollision
%probability
if prm.recoll_prob %adjust leaf optical properties as a function of clumping
    ratio_rho=rho/omega;
    ratio_tau=tau/omega;
    rho_=rho;
    tau_=tau;
    omega_=omega;
    for iv = canopy.nbot:canopy.ntop
        A_sh(iv)=(1-omega_)/(1-((1-canopy.clumpfac(iv))*omega_));
        omega_sh(iv)=1-A_sh(iv);
        rho(iv)=ratio_rho*omega_sh(iv);
        tau(iv)=ratio_tau*omega_sh(iv);
        omega(iv)=omega_sh(iv);
        omega_initial(iv)=omega_;

    end
else %create vertically constant profile for optical properties
    rho_=rho;
    tau_=tau;
    omega_=omega;
    for iv = canopy.nbot:canopy.ntop
        rho(iv)=rho_;
        tau(iv)=tau_;
        omega(iv)=omega_;
        omega_initial(iv)=omega_;
    end
end




% Cumulative leaf area index (from canopy top) at mid-layer

for iv = canopy.ntop: -1: canopy.nbot
    if (iv == canopy.ntop)
        canopy.sumlai(iv) = 0.5 * canopy.dlai(iv);
    else
        canopy.sumlai(iv) = canopy.sumlai(iv+1) + canopy.dlai(iv);
    end
end




% --- Atmospheric solar radiation. Solar radiation is given as a unit of visible radiation
% and a unit of near-infrared radiation, both split into direct and diffuse components.

for p = 1:params.npts
    atmos.solar_zenith(p) = sunang.theta_rad(p); %30 * (pi / 180);     % Solar zentih angle (radians)
    if atmos.solar_zenith(p)>1.544616, atmos.solar_zenith(p)=1.544616;end %this to avoid imaginary numbers when sun has not yet risen
    %sunang.theta_rad=atmos.solar_zenith(p);  %theta is zenith angle
    atmos.swskyb(p) = rad.inbeam(p);            % Direct beam solar radiation for X waveband (W/m2)
    atmos.swskyd(p) = rad.indiffuse(p);            % Diffuse solar radiation for X waveband (W/m2)
end

for p = 1:params.npts  %MB: needs to loop over all layers **************

    for iv = canopy.nbot:canopy.ntop

        G(iv)=0;
        for i=1:9
            x_=acos(cot(atmos.solar_zenith(p))*cot(theta_l(i)));
            if atmos.solar_zenith(p) <= pi/2-theta_l(i), S=cos(atmos.solar_zenith(p))*cos(theta_l(i));
            else S=cos(atmos.solar_zenith(p))*cos(theta_l(i))*(1+(2*(tan(x_)-x_)/pi));
            end
            G(iv)=G(iv)+g(leafanglecode(iv),i)*S;
        end
        G_wood(iv)=0;
        for i=1:9
            x_=acos(cot(atmos.solar_zenith(p))*cot(theta_l(i)));
            if atmos.solar_zenith(p) <= pi/2-theta_l(i), S=cos(atmos.solar_zenith(p))*cos(theta_l(i));
            else S=cos(atmos.solar_zenith(p))*cos(theta_l(i))*(1+(2*(tan(x_)-x_)/pi));
            end
            G_wood(iv)=G_wood(iv)+g_wood(leafanglecode_wood(iv),i)*S;
        end


        Kb(p,iv) = G(iv) / cos(atmos.solar_zenith(p));
        Kb_wood(p,iv) = G_wood(iv) / cos(atmos.solar_zenith(p));

    end


    %MB: cumulative effective LAI
    C_eff_LAI_times_K_b(p,canopy.ntop)=canopy.dlai(canopy.ntop)*canopy.clumpfac(canopy.ntop)*Kb(p,canopy.ntop);
    for iv = canopy.ntop-1: -1: canopy.nbot
        C_eff_LAI_times_K_b(p,iv)=C_eff_LAI_times_K_b(p,iv+1)+canopy.dlai(iv)*canopy.clumpfac(iv)*Kb(p,iv);
    end
    %MB: cumulative effective element area index (EAI)
    C_eff_EAI_times_K_b(p,canopy.ntop)=(canopy.dlai(canopy.ntop)*canopy.clumpfac(canopy.ntop)*Kb(p,canopy.ntop))+...
        (canopy.dwai(canopy.ntop)*canopy.clumpfac_wood(canopy.ntop)*Kb_wood(p,canopy.ntop));
    for iv = canopy.ntop-1: -1: canopy.nbot
        C_eff_EAI_times_K_b(p,iv)=C_eff_EAI_times_K_b(p,iv+1)+(canopy.dlai(iv)*canopy.clumpfac(iv)*Kb(p,iv))+...
            (canopy.dwai(iv)*canopy.clumpfac_wood(iv)*Kb_wood(p,iv));
    end
    %MB: cumulative effective wood area index (WAI)
    C_eff_WAI_times_K_b(p,canopy.ntop)=(canopy.dwai(canopy.ntop)*canopy.clumpfac_wood(canopy.ntop)*Kb_wood(p,canopy.ntop));
    for iv = canopy.ntop-1: -1: canopy.nbot
        C_eff_WAI_times_K_b(p,iv)=C_eff_WAI_times_K_b(p,iv+1)+canopy.dwai(iv)*canopy.clumpfac_wood(iv)*Kb_wood(p,iv);
    end


    if rad.incoming(p) == 0
        fsun(canopy.nbot:canopy.ntop)=0;
        fsun_wood(canopy.nbot:canopy.ntop)=0;
    else


        fsun(canopy.ntop)=(1-exp(-C_eff_EAI_times_K_b(p,canopy.ntop)))/(canopy.dlai(canopy.ntop)*Kb(p,canopy.ntop)+canopy.dwai(canopy.ntop)*Kb_wood(p,canopy.ntop));

        for iv = canopy.ntop-1: -1: canopy.nbot

            fsun(iv)=(exp(-C_eff_EAI_times_K_b(p,iv+1))-exp(-C_eff_EAI_times_K_b(p,iv)))/(canopy.dlai(iv)*Kb(p,iv)+canopy.dwai(iv)*Kb_wood(p,iv));

        end

        fsun_wood=fsun; %assigning same fsun to wood, as per discussion with Gordon on Nov 28 2023

        %removing very small values of fsun, this is causing problems in
        %longwave absorption because of a division by LAI*fsun
        fsun(fsun<0.0001)=0;
        fsun_wood(fsun_wood<0.0001)=0;
    end
    flux.fracsun(p,:)=fsun;
    flux.fracsun_wood(p,:)=fsun_wood;
    % --- Sunlit and shaded portions of canopy


    % Sunlit and shaded fraction of leaf layer

    for iv = canopy.nbot:canopy.ntop

        flux.fracsha(p,iv) = 1 - flux.fracsun(p,iv);
        flux.fracsha_wood(p,iv) = 1 - flux.fracsun_wood(p,iv);
    end


    for iv = canopy.nbot:canopy.ntop

        tb(p,iv) = exp(-(Kb(p,iv) * canopy.dlai(iv) * canopy.clumpfac(iv)+Kb_wood(p,iv) * canopy.dwai(iv) * canopy.clumpfac_wood(iv)));
        frac_b_l(p,iv)= (G(iv) * canopy.dlai(iv) * canopy.clumpfac(iv)) / (G(iv) * canopy.dlai(iv) * canopy.clumpfac(iv) + G_wood(iv) * canopy.dwai(iv) * canopy.clumpfac_wood(iv));
        frac_b_w(p,iv)= (G_wood(iv) * canopy.dwai(iv) * canopy.clumpfac_wood(iv)) / (G(iv) * canopy.dlai(iv) * canopy.clumpfac(iv) + G_wood(iv) * canopy.dwai(iv) * canopy.clumpfac_wood(iv));

    end

    for iv = canopy.nbot:canopy.ntop
        td(p,iv) = 0;

        frac_d_l(p,iv)=0;
        frac_d_w(p,iv)=0;

        for j = 1:90 %below I am not using sun zenith , but theta j which is G function when light is coming from different directions
            j_rad=j* pi / 180;
            G=0;
            for i=1:9
                x_=acos(cot(j_rad)*cot(theta_l(i)));
                if j_rad <= pi/2-theta_l(i), S=cos(j_rad)*cos(theta_l(i));
                else S=cos(j_rad)*cos(theta_l(i))*(1+(2*(tan(x_)-x_)/pi));
                end
                G=G+g(leafanglecode(iv),i)*S;
            end
            G_wood=0;
            for i=1:9
                x_=acos(cot(j_rad)*cot(theta_l(i)));
                if j_rad <= pi/2-theta_l(i), S=cos(j_rad)*cos(theta_l(i));
                else S=cos(j_rad)*cos(theta_l(i))*(1+(2*(tan(x_)-x_)/pi));
                end
                G_wood=G_wood+g_wood(leafanglecode_wood(iv),i)*S;
            end

            td(p,iv) = td(p,iv) ...
                + exp(-(G / cos(j_rad) * canopy.dlai(iv) * canopy.clumpfac(iv)+...
                G_wood / cos(j_rad) * canopy.dwai(iv) * canopy.clumpfac_wood(iv))) * sin(j_rad) * cos(j_rad);
            frac_d_l(p,iv)=frac_d_l(p,iv) + (G * canopy.dlai(iv) * canopy.clumpfac(iv)) / (G * canopy.dlai(iv) * canopy.clumpfac(iv) + G_wood * canopy.dwai(iv) * canopy.clumpfac_wood(iv)) * sin(j_rad) * cos(j_rad);
            frac_d_w(p,iv)=frac_d_w(p,iv) + (G_wood * canopy.dwai(iv) * canopy.clumpfac_wood(iv)) / (G * canopy.dlai(iv) * canopy.clumpfac(iv) + G_wood * canopy.dwai(iv) * canopy.clumpfac_wood(iv)) * sin(j_rad) * cos(j_rad);

        end

        td(p,iv) = td(p,iv) * 2 * (pi / 180); %this last part is delta Z in 14.33

        frac_d_l(p,iv)=frac_d_l(p,iv) * 2 * (pi / 180);
        frac_d_w(p,iv)=1-frac_d_l(p,iv);
    end

    tbcum(p,:)=0;


    tbcum(p,canopy.ntop) = 1;
    for iv = canopy.ntop-1: -1: canopy.nsoi

        tbcum(p,iv) = tbcum(p,iv+1)* ...
            exp(-(Kb(p,iv+1) * canopy.dlai(iv+1)* canopy.clumpfac(iv+1)+Kb_wood(p,iv+1) * canopy.dwai(iv+1)* canopy.clumpfac_wood(iv+1)));
    end



    % --- Light profile through canopy

    switch light
        case 'Norman'
            [flux] = NormanRadiation_mb_canveg_with_wood_for_release (rho, tau, omega, omega_initial, rho_wood, td, frac_d_l, frac_d_w, tb, frac_b_l, frac_b_w, tbcum, params, canopy, atmos, flux,p,Kb(p,:));

    end

    rad.diffuse(p,:)=flux.diffuse(p,:); 
    rad.direct(p,:)=flux.direct(p,:);
    rad.diffuse_wood(p,:)=flux.diffuse_wood(p,:);
    rad.direct_wood(p,:)=flux.direct_wood(p,:);
    rad.sun(p,:)=flux.sun(p,:);
    rad.shade(p,:)=flux.shade(p,:);
    rad.sun_wood(p,:)=flux.sun_wood(p,:);
    rad.shade_wood(p,:)=flux.shade_wood(p,:);
    rad.up_flux(p,:)=flux.swup;
    rad.dn_flux(p,:)=flux.swdn;

    rad.albedo(p)=flux.albcan(p);
    rad.soil_abs(p)=flux.swsoi(p);

end


rad.Kb=Kb;
rad.Kb_wood=Kb_wood;

rad.C_eff_EAI_times_K_b=C_eff_EAI_times_K_b;

rad.sun_abs=flux.swleaf(:,2:end,params.sun);
rad.sun_abs(isinf(rad.sun_abs)|isnan(rad.sun_abs))=0;
rad.sh_abs=flux.swleaf(:,2:end,params.sha);  %in umol m-2 of leaf s-1, no (fev 2021): in W m-2
rad.sun_abs_stem=flux.swstem(:,2:end,params.sun);
rad.sun_abs_stem(isinf(rad.sun_abs_stem)|isnan(rad.sun_abs_stem))=0;
rad.sh_abs_stem=flux.swstem(:,2:end,params.sha);
rad.wood_abs=flux.wood_abs(:,2:end);
rad.wood_abs(isinf(rad.wood_abs)|isnan(rad.wood_abs))=0;



rad.prob_beam=flux.fracsun(:,2:end);
rad.prob_shade=flux.fracsha(:,2:end);
rad.fsun_wood=flux.fracsun_wood(:,2:end);
rad.fshade_wood=flux.fracsha_wood(:,2:end);

rad.beam_flux=atmos.swskyb(:) .* tbcum(:,:); %MB: this is the flux of direct light pentrating through layers

rad.td=td;
rad.frac_d_l=frac_d_l;
rad.frac_d_w=frac_d_w;

rad.tb=tb;
rad.frac_b_l=frac_b_l;
rad.frac_b_w=frac_b_w;


rad.tbcum=tbcum;


if (waveband == 'par')
    figure(213)
    clf
    hold on
end

if (waveband == 'nir')
    figure(214)
    clf
    hold on
end


plot(mean(rad.up_flux(:,1:prm.jtot),1),prm.sumlai,'LineWidth', 2);
ax=gca;
set(ax, 'ydir','reverse');
plot(mean(rad.dn_flux(:,1:prm.jtot),1),prm.sumlai,'LineWidth', 2);
plot(mean(rad.beam_flux(:,1:prm.jtot),1),prm.sumlai,'LineWidth', 2);
plot(mean(rad.sun_abs(:,1:prm.jtot),1),prm.sumlai,'LineWidth', 2);
plot(mean(rad.sh_abs(:,1:prm.jtot),1),prm.sumlai,'LineWidth', 2);



legend('up','down','beam', 'sun abs', 'sh abs','Location','best')
xlabel('Radiation Flux Density')
ylabel('Canopy Cumulative LAI')
title('Norman ', waveband)


end