%% Get runup
%% DO NOT RUN AS FULL SCRIPT - !!! RUN IN SECTIONS !!!
clear all
data_dir = 'C:\Users\alange\Desktop\DATA\DATA\Runup';
%% CoastalLens
disp('How many cameras are you processing?')
cam_num = str2double(string(inputdlg('How many cameras?')));
close all

answer = questdlg('Do you have a .mat R file?', 'R file', 'Yes', 'No', 'Yes');
switch answer
    case 'Yes'
        disp('Please select the R file you want to load in.')
        [temp_file, temp_file_path] = uigetfile(data_dir, 'R file');
        load(fullfile(temp_file_path, temp_file)); clear temp_file*
    case 'No'
        %% Get extrinsics
        for cc = 1:cam_num

            %% Get intrinsics
            disp('Load in Intrinsics.mat')
            [temp_file, temp_file_path] = uigetfile(data_dir, 'Intrinsics File');
            load(fullfile(temp_file_path, temp_file)); clear temp_file*
            aa = whos; id = contains({aa.name}, 'cameraParams');
            eval(['cameraParams = ' aa(id).name])

            R(cc).cameraParams = cameraParams; clear cameraParams* id aa

            %% Get image coordinates
            sprintf('Load in Camera %i frame with GCPs visible.', cc)
            [temp_file, temp_file_path] = uigetfile('*.tiff', 'Camera Frame with GCP');

            R(cc).I = undistortImage(imread(fullfile(temp_file_path, temp_file)), R(cc).cameraParams);
            image_fig = figure(1);clf
            [image_gcp] = select_image_gcp(R(cc).I, image_fig);

            [world_gcp] = select_target_gcp;

            %% Get extrinsics - CoastalLens
            % if cc == 1
            %     ids = [1:4 6];
            % else
            ids =[1:length(image_gcp)];
            % end
            [worldPose, inlineid] = estworldpose(image_gcp(ids,:),world_gcp(ids,:), ...
                R(cc).cameraParams.Intrinsics, 'MaxReprojectionError', 5);

            R(cc).image_gcp = image_gcp;
            R(cc).world_gcp = world_gcp;
            R(cc).worldPose = worldPose;
            R(cc).used_gcps = ids;
            hGCP = figure(cc);clf
            imshow(R(cc).I)
            hold on
            scatter(image_gcp(:,1), image_gcp(:,2), 100, 'r', 'filled')
            for ii = 1:length(image_gcp)
                text(image_gcp(ii,1)+25, image_gcp(ii,2)-25, ['GCP ' char(string(ii))], 'FontSize', 14, 'BackgroundColor', 'w')
            end % for ii = 1:length(image_gcp)
            iP = world2img(world_gcp(ids,:),pose2extr(worldPose),R(cc).cameraParams.Intrinsics);
            scatter(iP(:,1), iP(:,2), 50, 'y', 'LineWidth', 3)
            %% CIRN
            %
            % extrinsicsInitialGuess = [R(cc).worldPose.Translation(1), R(cc).worldPose.Translation(2), R(cc).worldPose.Translation(3)...
            %     deg2rad(20) deg2rad(85) deg2rad(-10)];
            % extrinsicsKnownsFlag= [0 0 0 0 0 0];  % [ x y z azimuth tilt swing]
            % gcpInd = ids%[1:7]
            %
            % xyz = world_gcp(gcpInd,[1 2 3]);  % N x 3 matrix with rows= N gcps, columns= x,y,z
            % UVd=reshape([image_gcp(gcpInd,[2 1])],2,length(xyz))'; % N x 2 matrix with rows=gcps, columns= U,V
            % %
            % [extrinsics extrinsicsError]= extrinsicsSolver(extrinsicsInitialGuess,extrinsicsKnownsFlag,intrinsics_CIRN,UVd,xyz);
            % % Display the results
            % disp(' ')
            % disp('Solved Extrinsics and NLinfit Error')
            % disp( [' x = ' num2str(extrinsics(1)) ' +- ' num2str(extrinsicsError(1))])
            % disp( [' y = ' num2str(extrinsics(2)) ' +- ' num2str(extrinsicsError(2))])
            % disp( [' z = ' num2str(extrinsics(3)) ' +- ' num2str(extrinsicsError(3))])
            % disp( [' azimuth = ' num2str(rad2deg(extrinsics(4))) ' +- ' num2str(rad2deg(extrinsicsError(4))) ' degrees'])
            % disp( [' tilt = ' num2str(rad2deg(extrinsics(5))) ' +- ' num2str(rad2deg(extrinsicsError(5))) ' degrees'])
            % disp( [' swing = ' num2str(rad2deg(extrinsics(6))) ' +- ' num2str(rad2deg(extrinsicsError(6))) ' degrees'])
            % %
            % % % % Section 7: Reproject GCPs into UVd Space
            % xyzCheck = world_gcp(:,[1 2 3]);  % N x 3 matrix with rows= N gcps, columns= x,y,z
            % % Transform xyz World Coordinates to Distorted Image Coordinates
            % [UVdReproj ] = xyz2DistUV(intrinsics_CIRN,extrinsics,xyzCheck);
            % %  Reshape UVdCheck so easier to interpret
            % UVdReproj = reshape(UVdReproj ,[],2);
            %
            % % % Load Specified Image and Plot Clicked and Reprojected UV GCP Coordinates
            % f1=figure(3);clf
            % imshow(R(cc).I)
            % hold on
            %
            % for k=1:length(image_gcp)
            %     % Clicked Values
            %     h1=plot(image_gcp(k,1), image_gcp(k,2),'ro','markersize',10,'linewidth',3);
            %     % text(image_gcp(k,1)+30,image_gcp(k,2),num2str(k),'color','r','fontweight','bold','fontsize',15)
            %
            %     % New Reprojected Values
            %     h2=plot(UVdReproj(k,1),UVdReproj(k,2),'yo','markersize',10,'linewidth',3);
            %     % text(UVdReproj(k,1)+30,UVdReproj(k,2),num2str(gcp(k).num),'color','y','fontweight','bold','fontsize',15)
            % end
            % legend([h1 h2],'Clicked UVd','Reprojected UVd')

        end
end
%% Plot rectified image - grid product
disp('Only input grid product.')
answer = questdlg('Do you have a .mat Products file?', 'Product file', 'Yes', 'No', 'Yes');
switch answer
    case 'Yes'
        disp('Please select the grid product file you want to load in.')
        [temp_file, temp_file_path] = uigetfile(data_dir, 'Product file');
        load(fullfile(temp_file_path, temp_file)); clear temp_file*

        if ~exist('Products', 'var')
            disp('Please create Products file.')
            [Products] = user_input_products(data_dir);
        end % if ~exist('Products', 'var')
    case 'No'
        [Products] = user_input_products(data_dir);
end % switch answer
clear answer

% for Marconi CACO-03
% Products.type = 'Grid';
% Products.frameRate = 1;
% Products.lat = 41.8926;
% Products.lon = -69.9632;
% Products.angle = 80;
% Products.xlim = [200 0];
% Products.ylim = [-200 100];
% Products.dx = 0.1;
% Products.dy = 0.1;
% Products.x = [];
% Products.y = [];
% Products.z = [];
Products.tide = 0;

for cc = 1:cam_num
    if isfield(Products, 'iP')
        Products = rmfield(Products, 'iP');
        Products = rmfield(Products, 'iP_u');
        Products = rmfield(Products, 'iP_v');
    end

    if isfield(Products, 'Irgb_2d')
        Products = rmfield(Products, 'Irgb_2d');
    end


    [xyz, localX, localY, Z, Eastings, Northings] = getCoords(Products);

    aa=xyz-[R(cc).worldPose.Translation(1) R(cc).worldPose.Translation(2) 0];

    Products.xyz = xyz;
    Products.localX = localX;
    Products.localY = localY;
    Products.Eastings = Eastings;
    Products.Northings = Northings;
    Products.localZ = Z;

    Products.iP = round(world2img(Products.xyz, pose2extr(R(cc).worldPose), R(cc).cameraParams.Intrinsics));
    iP_u = reshape(Products.iP(:,2), size(Products.localX,1), size(Products.localX,2));
    iP_v = reshape(Products.iP(:,1), size(Products.localX,1), size(Products.localX,2));
    iP_u(iP_u <= 0) = NaN; iP_u(iP_u >= size(R(cc).I,1)) = NaN;
    iP_v(iP_v <= 0) = NaN; iP_v(iP_v >= size(R(cc).I,2)) = NaN;
    iP_u(isnan(iP_v)) = NaN; iP_v(isnan(iP_u)) = NaN;
    Products.iP_u = iP_u;
    Products.iP_v = iP_v;

    [rows, cols, numChannels] = size(R(cc).I);
    Irgb_temp = repmat(uint8([0]), size(Products.localX,1)*size(Products.localX,2),numChannels);

    for i = 1:numChannels
        channel = R(cc).I(:,:,i);
        Irgb_temp(~isnan(iP_u),i) = channel(sub2ind([rows, cols], iP_u(~isnan(iP_u)), iP_v(~isnan(iP_u))));
    end
    Irgb_temp=reshape(Irgb_temp, size(Products.localX,1),size(Products.localX,2),3);

    Products.Irgb_2d = uint8(Irgb_temp);

    IrIndv(:,:,:,cc) = squeeze(Products.Irgb_2d);

end % for cc = 1 : 2 % Cam 1 or 2


if sum(IrIndv(:)) ~= 0 % some color values present

    [Ir] =cameraSeamBlend(IrIndv);
    figure(1);clf
    image(Products(1).localX(:), Products(1).localY(:), Ir)
    axis equal
    xlim([min(Products(1).xlim) max(Products(1).xlim)])
    ylim([min(Products(1).ylim) max(Products(1).ylim)])
    set(gca, 'FontSize', 16)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.96]);
    xlabel('Cross-shore Distance (m)')
    ylabel('Along-shore Distance (m)')
    %  saveas(gcf, fullfile(data_dir, 'Processed_data', strcat('ARGUS_', day_files(dd).name, '_Grid.png')))
end

%% Get timestack coordinates



    disp('Only input timestacks.')
    answer = questdlg('Do you have a .mat Products file?', 'Product file', 'Yes', 'No', 'Yes');
    switch answer
        case 'Yes'
            disp('Please select the grid product file you want to load in.')
            [temp_file, temp_file_path] = uigetfile(data_dir, 'Product file');
            load(fullfile(temp_file_path, temp_file)); clear temp_file*
            
            if ~exist('Timestacks_CAM', 'var')
                if ~exist('Products', 'var')
                disp('Please create Timestacks_CAM file.')
                [Timestacks_CAM] = user_input_products(data_dir);
                else
                    Timestacks_CAM = Products(find(ismember(string({Products.type}), 'xTransect')));
                end

            end % if ~exist('Products', 'var')
        case 'No'
            [Timestacks_CAM] = user_input_products(data_dir);
    end % switch answer
    clear answer
    [Timestacks_CAM.tide]=deal(0);

  % for Marconi - CACO-03
    % y_opt = [-20 -30 -40 -50 -70 -100 -150];

for cc = 1:cam_num
  
    for yy = 1:length(Timestacks_CAM)
        % for Marconi - CACO-03
    %     Timestacks_CAM(yy).type = 'xTransect';
    %     Timestacks_CAM(yy).frameRate = 1;
    %     Timestacks_CAM(yy).lat = 41.8926;
    %     Timestacks_CAM(yy).lon = -69.9632;
    %     Timestacks_CAM(yy).angle = 80;
    %     Timestacks_CAM(yy).xlim = [200 0];
    %     Timestacks_CAM(yy).ylim = [];
    %     Timestacks_CAM(yy).dx = 0.1;
    %     Timestacks_CAM(yy).dy = 0.1;
    %     Timestacks_CAM(yy).x = [];
    %     Timestacks_CAM(yy).y = [];
    %     Timestacks_CAM(yy).z = [];
    %     Timestacks_CAM(yy).tide = 0;
    % 
    %     Timestacks_CAM(yy).y = y_opt(yy);

        [xyz, localX, localY, Z, Eastings, Northings] = getCoords(Timestacks_CAM(yy));

        aa=xyz-[R(cc).worldPose.Translation(1) R(cc).worldPose.Translation(2) 0];

        Timestacks_CAM(yy).xyz = xyz;
        Timestacks_CAM(yy).localX = localX;
        Timestacks_CAM(yy).localY = localY;
        Timestacks_CAM(yy).Eastings = Eastings;
        Timestacks_CAM(yy).Northings = Northings;
        Timestacks_CAM(yy).localZ = Z;

        Timestacks_CAM(yy).iP = round(world2img(Timestacks_CAM(yy).xyz, pose2extr(R(cc).worldPose), R(cc).cameraParams.Intrinsics));

        iP_u = Timestacks_CAM(yy).iP(:,1);
        iP_v = Timestacks_CAM(yy).iP(:,2);
        iP_u(iP_u <= 0) = NaN; iP_u(iP_u >= size(R(cc).I,2)) = NaN;
        iP_v(iP_v <= 0) = NaN; iP_v(iP_v >= size(R(cc).I,1)) = NaN;

        Timestacks_CAM(yy).localX_nan = Timestacks_CAM(yy).localX;
        Timestacks_CAM(yy).localY_nan = Timestacks_CAM(yy).localY;
        Timestacks_CAM(yy).localX_nan(isnan(iP_v))=NaN;
        Timestacks_CAM(yy).localX_nan(isnan(iP_u))=NaN;
        Timestacks_CAM(yy).localY_nan(isnan(iP_v))=NaN;
        Timestacks_CAM(yy).localY_nan(isnan(iP_u))=NaN;
        iP_u(isnan(iP_v)) = []; iP_v(isnan(iP_v)) = [];
        iP_v(isnan(iP_u)) = []; iP_u(isnan(iP_u)) = [];
        Timestacks_CAM(yy).iP = [];
        Timestacks_CAM(yy).iP = [iP_u iP_v];

    end

    eval(['Timestacks_CAM' char(string(cc)) '= Timestacks_CAM;'])
       
end % for cc = 1:cam_num

%% Make .pix files for i2R system
for cc = 1:cam_num
    eval(['c' char(string(cc)) '_timestack=[];'])
    eval(['for yy = 1:length(Timestacks_CAM' char(string(cc)) '); c' char(string(cc)) '_timestack = [c' char(string(cc)) '_timestack; Timestacks_CAM' char(string(cc)) '(yy).iP]; end'])
    eval(['fid = fopen(''c' char(string(cc)) '_timestack.pix'', ''w'');fprintf(fid, ''%i %i\n'', c' char(string(cc)) '_timestack'')'])
end
%% Import timestack coordinates from I2R

% scp -P 5000 argus_user@166.184.200.178:/mnt/I2Rgus_Data/ImageProducts/*.ras.tiff .
% put raw timestacks from i2R into 'Timestacks_raw' folder
data_files = dir(fullfile(data_dir, 'Timestacks_raw')); data_files(strcmp({data_files.name},{'.'}))=[]; data_files(strcmp({data_files.name},{'..'}))=[];

for cc = 1:cam_num
    eval(['data_files_c' char(string(cc)) ' = data_files(contains({data_files.name}, {''c' char(string(cc)) '''}));'])
end

for ll = 1:length(data_files_c1)
    clearvars -except data_dir data_files* ll cam_num Timestacks_CAM* Products Ir
    for cc = 1:cam_num
        clear Timestacks_CAM
        eval(['Timestacks_CAM = Timestacks_CAM' char(string(cc)) ';'])
        eval(['data_files_c = data_files_c' char(string(cc)) ';'])

       
        y_opt = [Timestacks_CAM.y];

        clear ts
        ts = imread(fullfile(data_files_c(ll).folder, data_files_c(ll).name));

        x=[Timestacks_CAM.localX];x = x(:); x(isnan(x))=[];
        for yy = 1:length(Timestacks_CAM)
            x_id = length(find(~isnan(Timestacks_CAM(yy).localX_nan)==1));
            try
                Timestacks_CAM(yy).Irgb_2d_short = ts(:, 1:x_id,:);

                ts(:, 1:x_id,:)=[];
                %image(Timestacks_CAM(yy).localX_nan(~isnan(Timestacks_CAM(yy).localX_nan)), [1:600], Timestacks_CAM(yy).Irgb_2d_short)

                Timestacks_CAM(yy).Irgb_2d = uint8(zeros(size(Timestacks_CAM(yy).Irgb_2d_short,1), size(Timestacks_CAM(yy).localX,1),3));
                id = find(Timestacks_CAM(yy).localX_nan == Timestacks_CAM(yy).localX);
                Timestacks_CAM(yy).Irgb_2d(:, id, :) = Timestacks_CAM(yy).Irgb_2d_short;
            end
            %nexttile()
            %image(Timestacks_CAM(yy).localX, [1:600], Timestacks_CAM(yy).Irgb_2d)
        end % for yy

        eval(['Timestacks_CAM' char(string(cc)) '= Timestacks_CAM;'])
        
    end % for cc = 1:cam_num

    %% combining transects - assumes each Timestacks_CAM structure has the same timestacks in it
    figure(ll);clf
    for yy = 1:length(Timestacks_CAM1)
        
        Image(yy).localX = Timestacks_CAM1(yy).localX;
        Image(yy).localY = Timestacks_CAM1(yy).localY;

        Image(yy).Eastings = Timestacks_CAM1(yy).Eastings;
        Image(yy).Northings = Timestacks_CAM1(yy).Northings;
        try
            for cc = 1:cam_num
                eval(['aa(:,:,:,' char(string(cc)) ') = Timestacks_CAM' char(string(cc)) '(yy).Irgb_2d;'])
            end
        [Image(yy).Irgb_2d] = cameraSeamBlend(aa);

        nexttile()
        image(Image(yy).localX, [1:size(Image(yy).Irgb_2d,2)], Image(yy).Irgb_2d)
        aa = split(data_files_c1(ll).name, '.');
        if size(Image(yy).Irgb_2d,1) == 600
            frameRate = 1;
        elseif size(Image(yy).Irgb_2d,1) == 1200
            frameRate = 2;
        end
        name = strcat(aa(1), '_', aa(2), '_', aa(3), '_', aa(4), '_', aa(5), '_', aa(6), '_', aa(7), '_dx', string(Timestacks_CAM1(yy).dx*10), '_yy', string(y_opt(yy)), '_', string(frameRate), 'Hz.png');
        
        ts = Image(yy).Irgb_2d; ts = fliplr(im2gray(ts));
        
        imwrite(ts, fullfile(pwd, 'Timestacks', name))
        end
    end % for yy

end % for ll

%% plot location on grid product

figure(1);clf
image(Products(1).localX(:), Products(1).localY(:), Ir)
axis equal
xlim([min(Products(1).xlim) max(Products(1).xlim)])
ylim([min(Products(1).ylim) max(Products(1).ylim)])
set(gca, 'FontSize', 16)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.96]);
xlabel('Cross-shore Distance (m)')
ylabel('Along-shore Distance (m)')

hold on

y_opt = [Timestacks_CAM.y];
for yy = 1:length(y_opt)
    for cc = 1:cam_num
    eval(['image(Timestacks_CAM' char(string(cc)) '(yy).localX, y_opt(yy), uint8(max(Timestacks_CAM' char(string(cc)) '(yy).Irgb_2d)))'])
    end
end

%% Extract horizontal runup

% do image segmentation
% conda info --envs
% conda activate gym_gpu
% go to folder with seg_images_in_folder.py
% run 'python seg_images_in_folder.py'
%

images_rgb = imageDatastore(fullfile(data_dir, 'Timestacks'));
images_pred = imageDatastore(fullfile(data_dir, 'Timestacks', 'meta'));

        %% Extract runup lines from segmented images
        for viewId = 1:length(images_pred.Files)
            disp(viewId)
            clear r_edge r_line x_line I

            aa = split(images_pred.Files(viewId), '/');
            oname = aa{end}(1:end-12);
            aa = split(oname, '_');

            I = sum(readimage(images_pred, viewId),3);
            if any(sum(I,2)==0)
                x_line = NaN;
                r_line = NaN;
            else

                r_edge = edge(I);
                u = unique(I(:));
                max_diff = u(end)-0;
                for vv = 1:size(r_edge,1)
                    if ~any(diff(I(vv,:))==max_diff) % no black line going into runup line
                        bb = find(r_edge(vv,:)==1);
                        if ~isempty(bb)
                            r_line(vv)=bb(end)-1;
                        else
                            r_line(vv)=NaN;
                        end
                    else
                        r_line(vv)=NaN;
                    end
                end
                % remove any points with too big a jump
                r_line(find(abs(diff(r_line))> 25)+1) = NaN;
                % remove any points outside of 3 std of swash
                r_line(find(r_line > median(r_line, 'omitnan')+2*std(r_line, 'omitnan')))=NaN;
                r_line(find(r_line < median(r_line, 'omitnan')-2*std(r_line, 'omitnan')))=NaN;
                % remove any points not on the line
                for ii=1:length(r_line)
                    if ~isnan(r_line(ii)) && I(ii, r_line(ii)-5)== I(ii,r_line(ii)+5)
                        r_line(ii)=NaN;
                    end
                end
                % interpolate acorss NaNs
                x_line = [0:size(r_edge,1)-1]*str2double(aa{end}(1:end-1));
                r_line = interp1(x_line(~isnan(r_line)), r_line(~isnan(r_line)), x_line);

                for vv = 1:length(r_line)
                    if any(diff(I(vv,:))==max(I(vv,:))) % no black line going into runup line
                        r_line(vv)=NaN;
                    end
                end

            end % if completly black across

            % clf
            % imshow(readimage(images_pred, viewId))
            % hold on
            % plot(r_line, [1:size(r_edge,1)], 'Color', 'w', 'LineWidth', 3)

            save(fullfile(data_dir, 'Timestacks', 'meta', strcat(oname, '_runup_line.mat')), 'r_line', 'x_line')

            
            % get MOP number
            % aa = split(oname, '_');
            % y_id = str2double(aa{find(strcmp(aa, 'y'))+1}(1:end-1));
            % if strcmp(aa{2}, 'Torrey')
            %     mop_origin = 582;
            % end
            % mop_numbers(viewId) = mop_origin + y_id/100;

        end

        %% Plot on image
        for viewId = 1:length(images_pred.Files)
            aa = split(images_pred.Files(viewId), '\');
            oname = aa{end}(1:end-12);

            load(fullfile(data_dir, 'Timestacks', 'meta', strcat(oname, '_runup_line.mat')), 'r_line', 'x_line')
            figure(1);clf
            imshow(readimage(images_rgb, viewId))
            hold on
            plot(r_line, [1:length(r_line)], 'Color', 'w', 'LineWidth', 3)
            pause(0.5)
        end

        %% =============== Only keep relevant MOP files  ============================
        % 
        % if all(unique(round(mop_numbers)) < 1000)
        %     MOP_files(~contains({MOP_files.name}, strcat('00', string(unique(round(mop_numbers)))))) = [];
        % elseif all(unique(round(mop_numbers)) > 1000)
        %     MOP_files(~contains({MOP_files.name}, strcat('0', string(unique(round(mop_numbers)))))) = [];
        % end

        %% Extract DEM and translate horizontal runup to vertical runup
        for viewId = 1:length(images_pred.Files)

            %load(fullfile(MOP_files(contains({MOP_files.name}, string(round(mop_numbers(viewId))))).folder, MOP_files(contains({MOP_files.name}, string(round(mop_numbers(viewId))))).name), 'SM');
            aa = split(images_pred.Files(viewId), '\');
            oname = aa{end}(1:end-12);
            aa = split(oname, '_');
            load(fullfile(data_dir, 'Timestacks', 'meta', strcat(oname, '_runup_line.mat')), 'r_line', 'x_line');

            [~, id] = min(abs(datenum(datetime(strcat(aa{6}), 'InputFormat', 'dd-MMM-yyyy')) - [SM.FileDatenum]));
            runup(viewId).mop = mop_numbers(viewId);
            runup(viewId).survey_date = datetime(SM(id).FileDatenum, 'ConvertFrom', 'datenum');
            runup(viewId).video_date = datetime(strcat(aa{6}, aa{7}), 'InputFormat', 'dd-MMM-yyyyHH:mm:ss', 'TimeZone', 'UTC');
            [~,~,verified,~,~] = getNOAAtide(runup(viewId).video_date,runup(viewId).video_date+minutes(17), '8447435');
            runup(viewId).tide = mean(verified);
            %runup(viewId).tide=0

            runup(viewId).dem_x = SM(id).X1D;
            runup(viewId).dem_z = SM(id).Z1Dmedian;

            frameRate = aa{find(strcmp(aa, 'frameRate'))+1};
            runup(viewId).time = [runup(viewId).video_date:seconds(1/str2double(frameRate(1:end-2))):runup(viewId).video_date+seconds((length(r_line)-1)*1/str2double(frameRate(1:end-2)))];

            for ll = find(~isnan(r_line))
                X = InterX([runup(viewId).dem_x; runup(viewId).dem_z],[x_line(int32(r_line(ll))) x_line(int32(r_line(ll))); 0 5]);
                if ~isempty(X)
                    runup(viewId).x(ll) = X(1);
                    runup(viewId).z(ll) = X(2);
                else
                    runup(viewId).x(ll) = NaN;
                    runup(viewId).z(ll) = NaN;
                end
            end
            if isnan(r_line)
                runup(viewId).x = NaN;
                runup(viewId).z = NaN;
            end
            runup(viewId).x(find(runup(viewId).x==0))=NaN;
            runup(viewId).z(find(runup(viewId).z==0))=NaN;

            figure(1);clf
            plot(runup(viewId).dem_x, runup(viewId).dem_z)
            hold on
            if all(~isnan(runup(viewId).x))
                plot(runup(viewId).x, runup(viewId).z, '.')
            else
                plot(x_line(int32(r_line)), runup(viewId).tide, '.')
            end
            yline(runup(viewId).tide)
            title(sprintf('MOP = %.2f',  runup(viewId).mop))
            %pause(1)
        end
        [~,id]=sort([runup.mop]);
        runup = runup(id);


        oname = [day_files(dd).name '_' flights(ff).name];
        save(fullfile(data_dir, 'Processed_data', [oname '_Runup']),'runup', '-v7.3')




