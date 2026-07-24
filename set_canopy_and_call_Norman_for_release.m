%this program sets a canopy structure, leaf optical properties, and sun
%geometry and calls the Norman model


clear all;
close all;


%for waveband=["nir" "par"]
for waveband=["nir"]

    %plot_list=["MorganMonroe" "EMS"  "Pasoh" "SERC"];
    plot_list=["MorganMonroe"];
    for plot_name_list=plot_list

        prm.plot_name=plot_name_list;

        % set illumination conditions
        sun_zenith_angle=35;
        diffuse_fraction=0.2;
        %set paths to files for LAI WAI and clumping profiles

        path_omega='/Users/mbeland/Documents/Research/Professor/Norman model checks/omega_corrected_';
        path_lai='/Users/mbeland/Documents/Research/Professor/Norman model checks/xAI_';

        sunang.theta_rad=deg2rad(sun_zenith_angle);

        sunang.theta_deg=sunang.theta_rad*180/pi;
        quantum.indiffuse=1000*diffuse_fraction;
        quantum.inbeam=1000-quantum.indiffuse;
        quantum.incoming = quantum.inbeam+ quantum.indiffuse;
        quantum.fraction_diffuse=quantum.indiffuse./quantum.incoming;

        file_omega=strcat(path_omega,prm.plot_name,'.xlsx');

        file_lai_FV=strcat(path_lai,prm.plot_name,'_with_angle_codes.asc');

        fid_lai_FV=fopen(file_lai_FV);
        omega_table=readtable(file_omega);

        %read omega data
        omega_leaf=omega_table.omega_leaf_corrected_30zen;

        omega_wood=omega_table.omega_wood_constant;

        %set canopy strucutre
        layer=0;
        %for i=1:2, fgets(fid_lai_FV);end
        fgets(fid_lai_FV);
        %for layer=canopy.nbot:canopy.ntop

        while ~feof(fid_lai_FV)
            layer=layer+1;
            height=fscanf(fid_lai_FV, '%f', [1])+.15;
            %bidon=fscanf(fid_lai_FV, '%f', [3]);
            leaf_area_index=fscanf(fid_lai_FV, '%f', [1]);
            if height>1.35 && leaf_area_index<0.003, break, end
            %bidon=fscanf(fid_lai_FV, '%f', [1]);
            wood_area_index=fscanf(fid_lai_FV, '%f', [1]);
            %bidon=fscanf(fid_lai_FV, '%f', [3]);

            prm.height_profile(layer)=height;

            %generic G functions 1: erec	2:spheri	3:uni	4:plano
            leaf_ang_dist_num=fscanf(fid_lai_FV, '%d', [1]);
            switch leaf_ang_dist_num
                case 1
                    prm.leafangle(layer)="erectophile";
                case 2
                    prm.leafangle(layer)="spherical";
                case 3
                    prm.leafangle(layer)="uniform";
                case 4
                    prm.leafangle(layer)="planophile";
            end

            wood_ang_dist_num=fscanf(fid_lai_FV, '%d', [1]);
            switch wood_ang_dist_num
                case 1
                    prm.leafangle_wood(layer)="erectophile";
                case 2
                    prm.leafangle_wood(layer)="spherical";
                case 3
                    prm.leafangle_wood(layer)="uniform";
                case 4
                    prm.leafangle_wood(layer)="planophile";
            end

            if height >0.5
                prm.markov(layer)=omega_leaf(round(height));
                prm.markov_wood(layer)=omega_wood(round(height));
            else
                prm.markov(layer)=omega_leaf(1);
                prm.markov_wood(layer)=omega_wood(1);
            end

            LAI_profile(layer)=leaf_area_index;
            WAI_profile(layer)=wood_area_index;
        end

        prm.height_profile=shiftdim(prm.height_profile);
        prm.leafangle=shiftdim(prm.leafangle);
        prm.leafangle_wood=shiftdim(prm.leafangle_wood);
        prm.markov=shiftdim(prm.markov);
        prm.markov_wood=shiftdim(prm.markov_wood);
        LAI_profile=shiftdim(LAI_profile);
        WAI_profile=shiftdim(WAI_profile);


        figure (100)
        clf
        hold on
        plot(LAI_profile,prm.height_profile,'LineWidth',2)
        plot(WAI_profile,prm.height_profile,'LineWidth',2)
        xlabel('Area index')
        ylabel ('height')
        legend('LAI','WAI')

        figure (101)
        clf
        hold on
        plot(prm.markov,prm.height_profile,'LineWidth',2)
        plot(prm.markov_wood,prm.height_profile,'LineWidth',2)
        xlabel('omega')
        ylabel ('height')
        legend('leaf','wood')


        prm.LAI=sum(LAI_profile);   % LIDAR LAI for each layer

        prm.sumlai=prm.LAI-cumsum(LAI_profile);

        prm.veg_ht=max(prm.height_profile);  %35.3;  % vegetation canopy height, LIDAR

        prm.dff=LAI_profile;    % each layer has a different LAI increment
        prm.dff(prm.dff==0)=0.00001; %MB: this is to avoid NaN values for LAI=0 layers in Norman RT
        prm.dfw=WAI_profile;    % each layer has a different WAI increment
        prm.dfw(prm.dfw==0)=0.00001; %MB: this is to avoid NaN values for WAI=0 layers in Norman RT

        prm.nlayers=length(LAI_profile);   % number of layers from LIDAR..top level was 0

        prm.jtot=int32(prm.nlayers);    % number of layers..had to int them. int8 did not work for 145 layers
        prm.jktot=prm.jtot+1;    % number of layers plus 1 %MB: where is the soil? why 1 more layer? is the soil layer 1? In Bonan Norman layer 1 is soil


        switch prm.plot_name
            case "MorganMonroe"
                prm.par_reflect = .06;   %MB: from excel file %MB: used by Ross 1981 refl 6% trans 9%
                prm.par_trans = .07; %MB: from excel file
                prm.nir_reflect = 0.45;   %MB: from excel file % 0.36..based on field and Ocean Optics...wt with Planck Law. High leaf N and high reflected NIR
                prm.nir_trans = 0.45;     % MB: from excel file % UCD presentation shows NIR transmission is about the same as reflectance
            case "EMS"
                prm.par_reflect = .06;   %MB: from excel file %MB: used by Ross 1981 refl 6% trans 9%
                prm.par_trans = .07; %MB: from excel file
                prm.nir_reflect = 0.45;   %MB: from excel file % 0.36..based on field and Ocean Optics...wt with Planck Law. High leaf N and high reflected NIR
                prm.nir_trans = 0.45;     % MB: from excel file % UCD presentation shows NIR transmission is about the same as reflectance
            case "Pasoh"
                %         prm.par_reflect = .11;   %MB: from excel file %MB: used by Ross 1981 refl 6% trans 9%
                %         prm.par_trans = .06; %MB: from excel file
                %         prm.nir_reflect = 0.46;   %MB: from excel file % 0.36..based on field and Ocean Optics...wt with Planck Law. High leaf N and high reflected NIR
                %         prm.nir_trans = 0.33;     % MB: from excel file % UCD presentation shows NIR transmission is about the same as reflectance
                prm.par_reflect = .06;   % the runs for Pasoh were redone in nov 2024 following strong indication from Canveg that leaf optics were wrong for this site, I used the same leaf optics as MMS and EMS
                prm.par_trans = .07; %MB:
                prm.nir_reflect = 0.45;
                prm.nir_trans = 0.45;     %

            case "SERC"
                prm.par_reflect = .06;   %MB: from excel file %MB: used by Ross 1981 refl 6% trans 9%
                prm.par_trans = .05; %MB: from excel file
                prm.nir_reflect = 0.45;   %MB: from excel file % 0.36..based on field and Ocean Optics...wt with Planck Law. High leaf N and high reflected NIR
                prm.nir_trans = 0.45;     % MB: from excel file % UCD presentation shows NIR transmission is about the same as reflectance
        end

        prm.par_absorbed = (1. - prm.par_reflect - prm.par_trans);
        prm.par_reflect_wood=0.21;
        prm.par_soil_refl = 0.122;    %MB: from excel file

        %         optical properties NIR wave band

        prm.nir_absorbed = (1. - prm.nir_reflect - prm.nir_trans);
        prm.nir_reflect_wood=0.49;
        prm.nir_soil_refl = 0.214; % MB: from excel file

        for i=1:length(LAI_profile)
            switch char(prm.leafangle(i))
                case 'planophile'
                    prm.leafanglecode(i)=1;

                case 'spherical'
                    prm.leafanglecode(i)=2;

                case 'erectophile'
                    prm.leafanglecode(i)=3;

                case 'plagiophile'
                    prm.leafanglecode(i)=4;

                case 'extremophile'
                    prm.leafanglecode(i)=5;

                case 'uniform'
                    prm.leafanglecode(i)=6;
            end
            switch char(prm.leafangle_wood(i))
                case 'planophile'
                    prm.leafanglecode_wood(i)=1;

                case 'spherical'
                    prm.leafanglecode_wood(i)=2;

                case 'erectophile'
                    prm.leafanglecode_wood(i)=3;

                case 'plagiophile'
                    prm.leafanglecode_wood(i)=4;

                case 'extremophile'
                    prm.leafanglecode_wood(i)=5;

                case 'uniform'
                    prm.leafanglecode_wood(i)=6;
            end

        end

        prm.recoll_prob=1;

        %calling Norman
        [quantum_wood]=fRadTranCanopy_Matrix_tridiagonal_with_wood_for_release(sunang,quantum,waveband,prm);


        figure(102)
        clf
        hold on
        plot(quantum_wood.prob_beam,prm.height_profile)

        ylabel('Height')
        xlabel ('fsun')
        title(plot_name_list, 'fsun')



        figure()
        clf
        hold on

        plot(quantum_wood.sh_abs_stem.*prm.dfw'.*(1-quantum_wood.fsun_wood)+quantum_wood.sun_abs_stem.*prm.dfw'.*quantum_wood.fsun_wood,prm.height_profile) %this converts to m-2 of ground area


        ylabel('Height')
        xlabel ('Absorption by wood, W/m^2 (unit ground)')

        title(plot_name_list,strcat(upper(waveband), ', wood'))


        figure()
        clf
        hold on
        plot(quantum_wood.sun_abs.*prm.dff'.*quantum_wood.prob_beam,prm.height_profile)

        ylabel('Height')
        xlabel ('Absorption by sunlit leaves, W/m^2 (unit ground)')

        title(plot_name_list,strcat(upper(waveband), ', sunlit leaves'))


        figure()
        clf
        hold on
        plot(quantum_wood.sh_abs.*prm.dff'.*quantum_wood.prob_shade,prm.height_profile)


        ylabel('Height')
        xlabel ('Absorption by shaded leaves, W/m^2 (unit ground)')

        title(plot_name_list,strcat(upper(waveband), ', shaded leaves'))

        clearvars -except waveband plot_list

    end
end