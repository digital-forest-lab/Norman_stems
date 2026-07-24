function [flux] = NormanRadiation_mb_canveg_with_wood_for_Gordon (rho, tau, omega, omega_initial, rho_wood, td, frac_d_l, frac_d_w, tb, frac_b_l, frac_b_w, tbcum, params, canopy, atmos, flux,p,Kb)

% Compute solar radiation transfer through canopy using Norman (1979)
%from Bonan 2019 'Climate change and
% terrestrial ecosystem modeling'

% -----------------------------------------------------------------------
% Input
% rho            ! Leaf reflectance
% tau            ! Leaf transmittance
% omega          ! Leaf scattering coefficient
% td             ! Exponential transmittance of diffuse radiation through a single leaf layer
% tb             ! Exponential transmittance of direct beam radiation through a single leaf layer
% tbcum          ! Cumulative exponential transmittance of direct beam onto a canopy layer
% params.numrad  ! Number of wavebands
% params.npts    ! Number of grid points to process
% params.sun     ! Index for sunlit leaf
% params.sha     ! Index for shaded leaf
% canopy.ntop    ! Index for top leaf layer
% canopy.nbot    ! Index for bottom leaf layer
% canopy.nsoi    ! First canopy layer is soil
% canopy.dlai    ! Layer leaf area index (m2/m2)
% atmos.swskyb   ! Atmospheric direct beam solar radiation (W/m2)
% atmos.swskyd   ! Atmospheric diffuse solar radiation (W/m2)
% flux.fracsun   ! Sunlit fraction of canopy layer
% flux.fracsha   ! Shaded fraction of canopy layer
% flux.albsoib   ! Direct beam albedo of ground (soil)
% flux.albsoid   ! Diffuse albedo of ground (soil)
%
% Output
% flux.swleaf    ! Leaf absorbed solar radiation (W/m2 leaf)
% flux.swveg     ! Absorbed solar radiation, vegetation (W/m2)
% flux.swvegsun  ! Absorbed solar radiation, sunlit canopy (W/m2)
% flux.swvegsha  ! Absorbed solar radiation, shaded canopy (W/m2)
% flux.swsoi     ! Absorbed solar radiation, ground (W/m2)
% flux.albcan    ! Albedo above canopy
% -----------------------------------------------------------------------

% --- Set up tridiagonal matrix

iv = canopy.nsoi;
swup(iv) = 0;
swdn(iv) = 0;

for iv = canopy.nbot:canopy.ntop
    swup(iv) = 0;
    swdn(iv) = 0;

end

% There are two equations for each canopy layer and the soil. The first
% equation is the upward flux and the second equation is the downward flux.

m = 0; % Initialize equation index for tridiagonal matrix

% Soil: upward flux

iv = canopy.nsoi;
m = m + 1;
a(m) = 0;
b(m) = 1;
c(m) = -flux.albsoid;
d(m) = atmos.swskyb(p) * tbcum(p,iv) * flux.albsoib;

A=td(p,iv+1)+((1-td(p,iv+1))*frac_d_l(p,iv+1)*tau(iv+1)+(1-td(p,iv+1))*frac_d_w(p,iv+1)*rho_wood/2); %here A is for iv+1

B=(1-td(p,iv+1))*frac_d_l(p,iv+1)*rho(iv+1)+(1-td(p,iv+1))*frac_d_w(p,iv+1)*rho_wood/2; %here B is for iv+1

m = m + 1;

a(m) =-(B-A^2/B); %this is -ai
b(m) = 1;

c(m) = -(A/B); %this is -bi

d(m) = atmos.swskyb(p) * tbcum(p,iv+1) * ((1 - tb(p,iv+1)) *frac_b_l(p,iv+1) * (tau(iv+1) - rho(iv+1) * A/B)+(1-tb(p,iv+1))*frac_b_w(p,iv+1)*rho_wood/2*(1-A/B)); %this is di
%MB: in the above, d(m) is di in book, a is first element in matrix line
    %(-ai), and c(m) is third element (-bi).

% Leaf layers, excluding top layer

for iv = canopy.nbot:canopy.ntop-1

    A=td(p,iv)+((1-td(p,iv))*frac_d_l(p,iv)*tau(iv)+(1-td(p,iv))*frac_d_w(p,iv)*rho_wood/2); %here A is for iv
    B=(1-td(p,iv))*frac_d_l(p,iv)*rho(iv)+(1-td(p,iv))*frac_d_w(p,iv)*rho_wood/2; %here B is for iv
    
    m = m + 1;
    
    a(m)= -(A/B); %this is -ei
    b(m) = 1;
    
    c(m)= -(B-A^2/B); %this is -fi
    
    d(m) = atmos.swskyb(p) * tbcum(p,iv) * ((1 - tb(p,iv)) *frac_b_l(p,iv) * (rho(iv) - tau(iv) * A/B)+(1-tb(p,iv))*frac_b_w(p,iv)*rho_wood/2*(1-A/B)); %this is ci
    %MB: in the above, d(m) is ci in book, a is first element in matrix line
    %(-ei), and c(m) is third element (-fi). Also, because this refers to ei
    %and fi, and not ei+1 and fi+1 as in the book, the indices for td are i
    %and not i+1 as in the book
    
    % Downward flux
    
    A=td(p,iv+1)+((1-td(p,iv+1))*frac_d_l(p,iv+1)*tau(iv+1)+(1-td(p,iv+1))*frac_d_w(p,iv+1)*rho_wood/2); %here A is for iv+1
    B=(1-td(p,iv+1))*frac_d_l(p,iv+1)*rho(iv+1)+(1-td(p,iv+1))*frac_d_w(p,iv+1)*rho_wood/2; %here B is for iv+1
    
    m = m + 1;
    
    a(m) =-(B-A^2/B); %this is -ai
    b(m) = 1;
    
    c(m) = -(A/B); %this is -bi
    
    d(m) = atmos.swskyb(p) * tbcum(p,iv+1) * ((1 - tb(p,iv+1)) *frac_b_l(p,iv+1) * (tau(iv+1) - rho(iv+1) * A/B)+(1-tb(p,iv+1))*frac_b_w(p,iv+1)*rho_wood/2*(1-A/B)); %this is di
    %MB: in the above, d(m) is di in book, a is first element in matrix line
    %(-ai), and c(m) is third element (-bi).

    %MB: in the above, d is di in book
   
end

% Top canopy layer: upward flux

iv = canopy.ntop;

A=td(p,iv)+((1-td(p,iv))*frac_d_l(p,iv)*tau(iv)+(1-td(p,iv))*frac_d_w(p,iv)*rho_wood/2); %here A is for iv
B=(1-td(p,iv))*frac_d_l(p,iv)*rho(iv)+(1-td(p,iv))*frac_d_w(p,iv)*rho_wood/2; %here B is for iv

m = m + 1;

a(m)= -(A/B); %this is -ei
b(m) = 1;

c(m)= -(B-A^2/B); %this is -fi

d(m) = atmos.swskyb(p) * tbcum(p,iv) * ((1 - tb(p,iv)) *frac_b_l(p,iv) * (rho(iv)  - tau(iv) * A/B)+(1-tb(p,iv))*frac_b_w(p,iv)*rho_wood/2*(1-A/B)); %this is ci

% Top canopy layer: downward flux

m = m + 1;
a(m) = 0;
b(m) = 1;
c(m) = 0;
d(m) = atmos.swskyd(p);

% --- Solve tridiagonal equations for fluxes

[u] = tridiagonal_solver (a, b, c, d, m);

m = 0;

% Soil fluxes

iv = canopy.nsoi;
m = m + 1;
swup(iv) = u(m);
m = m + 1;
swdn(iv) = u(m);

% Leaf layer fluxes

for iv = canopy.nbot:canopy.ntop
    m = m + 1;
    swup(iv) = u(m);
    m = m + 1;
    swdn(iv) = u(m);
end

% --- Compute flux densities

% Absorbed direct beam and diffuse for ground (soil)

iv = canopy.nsoi;
direct = atmos.swskyb(p) * tbcum(p,iv) * (1 - flux.albsoib);
diffuse = swdn(iv) * (1 - flux.albsoid);
flux.swsoi(p) = direct + diffuse;

% Absorbed direct beam and diffuse for each leaf layer and sum
% for all leaf layers

flux.swveg(p) = 0;
flux.swvegsun(p) = 0;
flux.swvegsha(p) = 0;
flux.swwood(p) = 0;

for iv = canopy.nbot:canopy.ntop
  
    % Per unit ground area (W/m2 ground)

    diffuse = (swdn(iv) + swup(iv-1)) * (1 - td(p,iv)) * frac_d_l(p,iv) * (1 - omega(iv)); %eq 14.47
    diffuse = diffuse + atmos.swskyb(p) * tbcum(p,iv) * (1 - tb(p,iv)) * frac_b_l(p,iv) * (omega_initial(iv)-omega(iv)); %I'm affecting here the part of the beam absorption caused by recollision probability to diffuse light, because is is no longer beam radiation and does not contribute to direct light on sunlit leaves
    flux.diffuse(p,iv)=diffuse;
    
    diffuse_wood = (swdn(iv) + swup(iv-1)) * (1 - td(p,iv)) * frac_d_w(p,iv) * (1 - rho_wood); 
    flux.diffuse_wood(p,iv)=diffuse_wood;
    
    direct = atmos.swskyb(p) * tbcum(p,iv) * (1 - tb(p,iv)) * frac_b_l(p,iv) * (1 - omega_initial(iv)); %eq 14.48 Icb
    
    flux.direct(p,iv)=direct;
    
    direct_wood = atmos.swskyb(p) * tbcum(p,iv) * (1 - tb(p,iv)) * frac_b_w(p,iv) * (1 - rho_wood);
    flux.direct_wood(p,iv)=direct_wood;
    
    % Absorbed solar radiation for shaded and sunlit portions of leaf layer
    % per unit ground area (W/m2 ground)
    
    sun = diffuse * flux.fracsun(p,iv) + direct; %fracsun is fsun
    sun_wood = diffuse_wood * flux.fracsun_wood(p,iv) + direct_wood; 
    flux.sun(p,iv)=sun;
    flux.sun_wood(p,iv)=sun_wood;
    shade = diffuse * flux.fracsha(p,iv); %differs from 14.52
    shade_wood = diffuse_wood * flux.fracsha_wood(p,iv); 
    flux.shade(p,iv)=shade;
    flux.shade_wood(p,iv)=shade_wood;
    % Convert to per unit sunlit and shaded leaf area (W/m2 leaf)
    
    flux.swleaf(p,iv,params.sha) = shade / (flux.fracsha(p,iv) * canopy.dlai(iv)); %MB: 6 march 2021: this fits eq 14.52
    flux.swstem(p,iv,params.sha) = shade_wood / (flux.fracsha_wood(p,iv) * canopy.dwai(iv)); 
    
    flux.swleaf(p,iv,params.sun)= sun/ (flux.fracsun(p,iv) * canopy.dlai(iv)); 
    flux.swstem(p,iv,params.sun)= sun_wood/ (flux.fracsun_wood(p,iv) * canopy.dwai(iv)); %this is now per m2 wood
    
    flux.wood_abs(p,iv)=(sun_wood+shade_wood) / canopy.dwai(iv); %******this is the sum of absorbed light by shaded and sunlit wood converted to per m2 wood units
    flux.wood_abs_check(p,iv)=(direct_wood + diffuse_wood) / canopy.dwai(iv); %should be same as above
    
    % Sum fluxes over all leaf layers
    flux.swveg(p) = flux.swveg(p) + (direct + diffuse); %eq 14.49
    flux.swvegsun(p) = flux.swvegsun(p) + sun;
    flux.swvegsha(p) = flux.swvegsha(p) + shade;
    flux.swwood(p) = flux.swwood(p) + (direct_wood + diffuse_wood); %this is vertically integrated
end

% --- Albedo

incoming = atmos.swskyb(p) + atmos.swskyd(p);
reflected = swup(canopy.ntop); 
if (incoming > 0)
    flux.albcan(p) = reflected / incoming;
else
    flux.albcan(p) = 0;
end

% --- Conservation check

% Total radiation balance: absorbed = incoming - outgoing

suminc = atmos.swskyb(p) + atmos.swskyd(p);
sumref = flux.albcan(p) * (atmos.swskyb(p) + atmos.swskyd(p));
sumabs = suminc - sumref;

err = sumabs - (flux.swveg(p) + flux.swsoi(p) + flux.swwood(p));

if (abs(err) > 1e-03)
    fprintf('time step: %d\n',p)
    fprintf('err = %15.5f\n',err)
    fprintf('sumabs = %15.5f\n',sumabs)
    fprintf('swveg = %15.5f\n',flux.swveg(p))
    fprintf('swsoi = %15.5f\n',flux.swsoi(p))
    fprintf('swwood = %15.5f\n',flux.swwood(p))
    warning ('NormanRadiation: Total solar conservation error')
end

% Sunlit and shaded absorption

err = (flux.swvegsun(p) + flux.swvegsha(p)) - flux.swveg(p);
if (abs(err) > 1e-03)
    fprintf('err = %15.5f\n',err)
    fprintf('swveg = %15.5f\n',flux.swveg(p))
    fprintf('swvegsun = %15.5f\n',flux.swvegsun(p))
    fprintf('swvegsha = %15.5f\n',flux.swvegsha(p))
    warning ('NormanRadiation: Sunlit/shade solar conservation error')
end

%MB: adding variables output needed by Canveg

flux.swup=swup; %MB: revetred back because this was resulting in layer 1 being 0 in Canveg
flux.swdn=swdn;

end