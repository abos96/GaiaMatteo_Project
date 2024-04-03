classdef wacomdata
% Reads the tablet's data file and imports (corrected) data into Matlab
    properties
        time = [];
        orientation = 'l';  % landscape
        resolution = 166;
        fc = 200;
        position = [];
        Pm = [44800 29600]/2;   % Digitizer dimensions (midpoints)
        pressure = [];
        nsamples = 0;
    end

    methods
        function wd = read(wd,file,orient)
            data = readtable(file);

            t = datetime(char(data.Timestamp),InputFormat='uuuu-MM-dd''T''HH:mm:ss.SSSXXX',TimeZone='local');
            wd.time = seconds(t-t(1));

            wd.orientation = orient;

            switch wd.orientation
                case 'l'
                    X = data.PointX;
                    Y = -data.PointY;

                    wd.Pm(2) = -wd.Pm(2);
                case 'p'
                    X = data.PointY;
                    Y = data.PointX;

                    wd.Pm = [wd.Pm(2) wd.Pm(1)];
            end
            if Y(1)<(Y(1)+Y(end))/2
                X = -X;
                Y = -Y;

                wd.Pm = -[wd.Pm(1) wd.Pm(2)];
            end
            wd.position = [X Y]/wd.resolution;  % From px to mm
            wd.Pm = wd.Pm/wd.resolution;

            wd.pressure = data.Pressure;

            wd.nsamples = size(data,1);

            [wd.time,wd.position,wd.pressure] = data_correction1(wd); % Matteo - Gaia

            % [wd.time,wd.position,wd.pressure] = data_correction2(wd);   % Alessio
        end

        function [T,POS,PRES] = data_correction1(wd)
            Time = 0:1/wd.fc:wd.time(end);  % Corrected time array

            time1 = [0;wd.time(1:end-1)];
            TIME1 = wd.time-time1;  % Difference between raw-time and translated raw-time (with error in data capture)
            t = unique(wd.time,'stable');   % Delete repeated values in time
            time2 = [0;t(1:end-1)];
            TIME2 = t-time2;    % Difference between time and translated time (after deleting repeated values)

            % Find moments when there are (before) "repeated values" (idx1) and (afterwards) "holes" (idx2) during the data capture
            % We don't consider TIME(1) because it is 0, and this is not an acquisition error
            idx1 = find(round(TIME1(2:end),3)==0);
            idx2 = find(round(TIME2(2:end),3)~=1/wd.fc);

            wd.position(idx1,:) = [];   % Delete repeated values in position
            wd.pressure(idx1) = []; % Delete repeated values in pressure

            % First values are correct
            pos1 = wd.position(1:idx2(1),:);
            pres1 = wd.pressure(1:idx2(1));

            % Corrected array
            ps = [];
            prs = [];

            for i=1:length(idx2)-1
                n = round((t(idx2(i)+1)-t(idx2(i)))*wd.fc); % Number of replacement blocks
                N = idx2(i)+1:idx2(i+1);    % Indx of shifted values

                % Replacement blocks
                pos2{i} = nan(n,2);
                pres2{i} = zeros(n,1);

                % Shifted (correct) values
                pos3{i} = wd.position(N,:);
                pres3{i} = wd.pressure(N);

                ps = [ps;pos2{i};pos3{i}];
                prs = [prs;pres2{i};pres3{i}];
            end

            % Full corrected array
            ps = [pos1;ps;wd.position(idx2(end)+1:end,:)];
            prs = [pres1;prs;wd.pressure(idx2(end)+1:end)];

            % Interpolation of array
            tm = linspace(0,wd.time(end),length(ps));   % Time array with length of corrected arrays
            pos_int.x = interp1(tm,ps(:,1),Time,'nearest');
            pos_int.y = interp1(tm,ps(:,2),Time,'nearest');
            pres_int = interp1(tm,prs,Time,'nearest');

            % Eliminate additional jumps in position data
            x = pos_int.x-[0,pos_int.x(1:end-1)];

            for i=3:length(x)   % First and second values are not to be taken into account
                if x(i)-x(i-1)<-1 || x(i)-x(i-1)>1
                    pos_int.x(i) = NaN;
                    pos_int.y(i) = NaN;
                end
            end

            % Import the corrected data into the class
            T = Time';
            POS = [pos_int.x' pos_int.y'];
            PRES = pres_int';
        end

        % function [T,POS,PRES] = data_correction2(wd)
        %     % Create new time array
        %     Time = 0:1/wd.fc:wd.time(end);
        % 
        %     % Find not connected events
        %     idx = 0;
        % 
        %     for i=1:length(wd.time)-1
        %         if wd.time(i+1)-wd.time(i)>0.005001
        %             idx = idx+1;
        %             NaN_samples = round((wd.time(i+1)-wd.time(i))*wd.fc);
        %             NaN_index = floor(wd.time(i)*wd.fc);
        %             NaN_indeces{idx} = NaN_index:NaN_index+NaN_samples;
        %         end
        %     end
        % 
        %     % Interpolation
        %     % Find unique values in T
        %     [~,unique_indices] = unique(wd.time,'stable');
        % 
        %     % Delete not unique values
        %     t = wd.time(unique_indices);
        %     pos = wd.position(unique_indices,:);
        %     pres = wd.pressure(unique_indices);
        % 
        %     % Interpolation in the new time values
        %     pos_int.x = interp1(t,pos(:,1),Time,'pchip');
        %     pos_int.y = interp1(t,pos(:,2),Time,'pchip');
        %     pres_int = interp1(t,pres,Time,'pchip');
        % 
        %     % Add NaN and 0
        %     for i=1:length(NaN_indeces)
        %         pos_int.x(NaN_indeces{i}) = NaN;
        %         pos_int.y(NaN_indeces{i}) = NaN;
        %         pres_int(NaN_indeces{i}) = 0;
        %     end
        % 
        %     % Corrected data
        %     T = Time';
        %     POS = [pos_int.x' pos_int.y'];
        %     PRES = pres_int';
        % end
    end
end