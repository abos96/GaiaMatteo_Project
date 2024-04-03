classdef hw_text
    properties
        time = [];
        position = [];
        fs = [];
        SGorder = 4;
        SGwin = 11;
        pressure = [];
        incontact = [];
        inair = [];
        pres_mean = [];
        pres_std = [];
        vel = [];
        speed = [];
        speed_mean = [];
        speed_std = [];
        acc = [];
        acc_mean = [];
        acc_std = [];
        jer = [];
        jer_norm = [];
        jer_mean = [];
        jer_std = [];
        Rc = [];
        Rc_mean = [];
        Rc_std = [];
        Pc = [];
        Pa = [];
        % rw = [];
        % rows = [];
        % nRows = [];
        chunks = [];
        nChunks = [];
        words = [];
        nWords = [];
        theta = [];
        rho = [];
    end

    methods
        function txt = hw_text(wd,SGwin,SGorder)
            txt.time = wd.time;

            txt.position = wd.position;

            txt.fs = wd.fc;

            txt.SGorder = SGorder;

            txt.SGwin = SGwin;

            txt.pressure = wd.pressure;

            txt.incontact = find(txt.pressure>0);

            txt.inair = find(txt.pressure==0);

            [txt.pres_mean,txt.pres_std] = pressure_mean_std(txt);  % Calculate mean and std of pressure

            [txt.vel,txt.speed,txt.speed_mean,txt.speed_std] = velocity(txt);   % Calculate velocity, speed, mean and std of speed (in contact and in air)

            [txt.acc,txt.acc_mean,txt.acc_std] = acceleration(txt); % Calculate acceleration, mean and std of acceleration (in contact and in air)

            [txt.jer,txt.jer_norm,txt.jer_mean,txt.jer_std] = jerk(txt);    % Calculate jerk, jerk-norm, mean and std of jerk (in contact and in air)

            [txt.Rc,txt.Rc_mean,txt.Rc_std] = curvrad(txt); % Calculate curvature radius, its mean and std

            [txt.Pc,txt.Pa] = time_percentage(txt); % Calculate time percentage (in contact and in air)

            % [txt.rw,txt.rows,txt.nRows] = txt_segment(txt); % Segment text into rows

            [txt.chunks,txt.nChunks] = chunks_txt_segment(txt); % Segment text into chunks

            [txt.words,txt.nWords] = words_txt_segment(txt);    % Segment text into words

            [txt.theta,txt.rho] = words_polar_coord(txt,wd);   % Words distribution in polar coordinates
        end

        function plot_trajectory(txt)   % Trajectories plot (in air and in contact) and histogram of rows
            a5 = [210 148]*2.8346438836889; % Paper size A5

            fig = gcf;
            % fig.Position = [100 100 a5(1)+200 a5(2)];   % Change figure size
            fig.Position = [100 100 a5];

            % h1 = axes('outerpos',[0 0 0.15 1]);
            % h2 = axes('outerpos',[0.1 0 0.9 1]);

            % axes(h2)    % Plot text separator lines
            % for i=1:length(txt.rw.im)
            %     line([min(txt.position(txt.incontact,1)) max(txt.position(txt.incontact,1))],[txt.rw.HVyrange(txt.rw.im(i)) txt.rw.HVyrange(txt.rw.im(i))],'col','k','lines','-')
            % end

            % for i=1:txt.nRows   % Plot boxes of different lines
            %     line([txt.rows{i}.xmin txt.rows{i}.xmax],[txt.rows{i}.ymin txt.rows{i}.ymin],'col','k','lines','--')
            %     line([txt.rows{i}.xmin txt.rows{i}.xmax],[txt.rows{i}.ymax txt.rows{i}.ymax],'col','k','lines','--')
            % 
            %     line([txt.rows{i}.xmin txt.rows{i}.xmin],[txt.rows{i}.ymin txt.rows{i}.ymax],'col','k','lines','--')
            %     line([txt.rows{i}.xmax txt.rows{i}.xmax],[txt.rows{i}.ymin txt.rows{i}.ymax],'col','k','lines','--')
            % end

            % Trajectories plot
            % h(1) = line(txt.position(txt.incontact,1),txt.position(txt.incontact,2),'marker','.','col','b','lines','none');
            % h(2) = line(txt.position(txt.inair,1),txt.position(txt.inair,2),'marker','.','col','r','lines','none');
            % legend(h,'On-surface movement','In-air movement')

            line(txt.position(txt.incontact,1),txt.position(txt.incontact,2),'marker','.','col','b','lines','none');
            line(txt.position(txt.inair,1),txt.position(txt.inair,2),'marker','.','col','r','lines','none');
            legend('On-surface movement','In-air movement')

            axis equal
            grid on

            % yrange = get(h2,'ylim');    % Range of y-values
            % 
            % axes(h1)    % Plot of histogram
            % line(txt.rw.HVfilt,txt.rw.HVyrange)
            % 
            % for i=1:length(txt.rw.im)   % Plot of markers and lines from histogram
            %     line(-txt.rw.m(i),txt.rw.HVyrange(txt.rw.im(i)),'col','k','marker','o')
            %     line([-txt.rw.m(i) max(txt.rw.HVfilt)],[txt.rw.HVyrange(txt.rw.im(i)) txt.rw.HVyrange(txt.rw.im(i))],'col','k','lines','-.')
            % end
            % 
            % ylim(yrange)
        end

        function plot_trajectory_incontact(txt) % Trajectory plot (only in contact)
            a5 = [210 148]*2.8346438836889;

            fig = gcf;
            fig.Position = [100 100 a5];

            line(txt.position(txt.incontact,1),txt.position(txt.incontact,2),'marker','.','col','b','lines','none');
        end

        function [pres_mean,pres_std] = pressure_mean_std(txt)
            pres_mean = mean(txt.pressure(txt.incontact));
            pres_std = std(txt.pressure(txt.incontact));
        end

        function plot_pressure(txt) % Plot pressure and mean pressure
            line(txt.time,txt.pressure,'col','b')
            line([txt.time(1) txt.time(end)],[txt.pres_mean txt.pres_mean],'col','r')
            xlabel('Time [s]')
            % ylabel('P [MPa]')
            ylabel('P levels')
            title('Pressure')
            legend('Pressure','Mean pressure')
        end

        function plot_pressure_distribution(txt)
            pressure_no_outliers = rmoutliers(txt.pressure);    % Outliers removal

            boxplot(pressure_no_outliers)
            title('Pressure distribution plot')
            xticklabels('Pressure')
        end

        function [vel,speed,speed_mean,speed_std] = velocity(txt)
            vx = ((-txt.fs)^1)*savitzkyGolayFilt(txt.position(:,1),txt.SGorder,1,txt.SGwin);
            vy = ((-txt.fs)^1)*savitzkyGolayFilt(txt.position(:,2),txt.SGorder,1,txt.SGwin);

            vel.incontact = nan(length(txt.time),2);
            vel.incontact(txt.incontact,:) = [vx(txt.incontact) vy(txt.incontact)];

            vel.inair = nan(length(txt.time),2);
            vel.inair(txt.inair,:) = [vx(txt.inair) vy(txt.inair)];

            speed.incontact = nan(length(txt.time),1);
            speed.incontact(txt.incontact) = sqrt(vx(txt.incontact).^2+vy(txt.incontact).^2);

            speed.inair = nan(length(txt.time),1);
            speed.inair(txt.inair) = sqrt(vx(txt.inair).^2+vy(txt.inair).^2);

            speed_mean.incontact = mean(speed.incontact,"omitnan");
            speed_std.incontact = std(speed.incontact,"omitnan");

            speed_mean.inair = mean(speed.inair,"omitnan");
            speed_std.inair = std(speed.inair,"omitnan");
        end

        function plot_speed_incontact(txt)
            line(txt.time,txt.speed.incontact)
            xlabel('Time [s]')
            ylabel('V [mm/s]')
            title('Speed in contact')
        end

        function plot_speed_inair(txt)
            line(txt.time,txt.speed.inair)
            xlabel('Time [s]')
            ylabel('V [mm/s]')
            title('Speed in air')
        end

        function plot_speed_distribution(txt)
            speed_no_outliers = rmoutliers([txt.speed.incontact txt.speed.inair]);

            boxplot(speed_no_outliers)
            title('Speed distribution plot')
            xticklabels({'Speed in contact','Speed in air'})
        end

        function [acc,acc_mean,acc_std] = acceleration(txt)
            ax = ((-txt.fs)^2)*savitzkyGolayFilt(txt.position(:,1),txt.SGorder,2,txt.SGwin);
            ay = ((-txt.fs)^2)*savitzkyGolayFilt(txt.position(:,2),txt.SGorder,2,txt.SGwin);

            acc.incontact = nan(length(txt.time),2);
            acc.incontact(txt.incontact,:) = [ax(txt.incontact) ay(txt.incontact)];

            acc.inair = nan(length(txt.time),2);
            acc.inair(txt.inair,:) = [ax(txt.inair) ay(txt.inair)];

            acc_mean.incontact = mean(acc.incontact,"omitnan");
            acc_std.incontact = std(acc.incontact,"omitnan");

            acc_mean.inair = mean(acc.inair,"omitnan");
            acc_std.inair = std(acc.inair,"omitnan");
        end

        function plot_acceleration_incontact(txt)
            subplot(2,1,1)
            line(txt.time,txt.acc.incontact(:,1),'col','b')
            ylabel('A_x [mm/s^2]')
            title('Acceleration in contact')
            subplot(2,1,2)
            line(txt.time,txt.acc.incontact(:,2),'col','b')
            xlabel('Time [s]')
            ylabel('A_y [mm/s^2]')
        end

        function plot_acceleration_inair(txt)
            subplot(2,1,1)
            line(txt.time,txt.acc.inair(:,1),'col','b')
            ylabel('A_x [mm/s^2]')
            title('Acceleration in air')
            subplot(2,1,2)
            line(txt.time,txt.acc.inair(:,2),'col','b')
            xlabel('Time [s]')
            ylabel('A_y [mm/s^2]')
        end

        function plot_acceleration_distribution(txt)
            acc_no_outliers = rmoutliers([txt.acc.incontact txt.acc.inair]);

            boxplot(acc_no_outliers)
            title('Acceleration distribution plot')
            xticklabels({'Ax in contact','Ay in contact','Ax in air','Ay in air'})
        end

        function [jer,jer_norm,jer_mean,jer_std] = jerk(txt)
            jx = ((-txt.fs)^3)*savitzkyGolayFilt(txt.position(:,1),txt.SGorder,3,txt.SGwin);
            jy = ((-txt.fs)^3)*savitzkyGolayFilt(txt.position(:,2),txt.SGorder,3,txt.SGwin);

            jer.incontact = nan(length(txt.time),2);
            jer.incontact(txt.incontact,:) = [jx(txt.incontact) jy(txt.incontact)];

            jer.inair = nan(length(txt.time),2);
            jer.inair(txt.inair,:) = [jx(txt.inair) jy(txt.inair)];

            jer_norm.incontact = nan(length(txt.time),1);
            jer_norm.incontact(txt.incontact) = sqrt(jx(txt.incontact).^2+jy(txt.incontact).^2);

            jer_norm.inair = nan(length(txt.time),1);
            jer_norm.inair(txt.inair) = sqrt(jx(txt.inair).^2+jy(txt.inair).^2);

            jer_mean.incontact = mean(jer_norm.incontact,"omitnan");
            jer_std.incontact = std(jer_norm.incontact,"omitnan");

            jer_mean.inair = mean(jer_norm.inair,"omitnan");
            jer_std.inair = std(jer_norm.inair,"omitnan");
        end

        function plot_jerk_incontact(txt)
            line(txt.time,txt.jer_norm.incontact)
            xlabel('Time [s]')
            ylabel('J [mm/s^3]')
            title('Jerk in contact')
        end

        function plot_jerk_inair(txt)
            line(txt.time,txt.jer_norm.inair)
            xlabel('Time [s]')
            ylabel('J [mm/s^3]')
            title('Jerk in air')
        end

        function plot_jerk_distribution(txt)
            jerk_no_outliers = rmoutliers([txt.jer_norm.incontact txt.jer_norm.inair]);

            boxplot(jerk_no_outliers)
            title('Jerk distribution plot')
            xticklabels({'Jerk in contact','Jerk in air'})
        end

        function [Rc,Rc_mean,Rc_std] = curvrad(txt)
            % Formula from Wikipedia
            k = abs(txt.vel.incontact(:,1).*txt.acc.incontact(:,2)-txt.vel.incontact(:,2).*txt.acc.incontact(:,1))./((txt.vel.incontact(:,1).^2+txt.vel.incontact(:,2).^2).^(3/2));
            Rc = 1./k;

            Rc_mean = mean(Rc,"omitnan");
            Rc_std = std(Rc,"omitnan");
        end

        function plot_curvrad(txt)
            line(txt.time,txt.Rc)
            xlabel('Time [s]')
            ylabel('Rc [mm]')
            title('Curvature radius')
            set(gca,'YScale','log')
        end

        function plot_curvrad_distribution(txt)
            curvrad_no_outliers = rmoutliers(txt.Rc);

            boxplot(curvrad_no_outliers)
            title('Curvature radius distribution plot')
            xticklabels('Curvature radius')
        end

        function [Pc,Pa] = time_percentage(txt)
            Pc = length(txt.time(txt.incontact))/length(txt.time)*100;
            Pa = length(txt.time(txt.inair))/length(txt.time)*100;
        end

        % function [rw,rows,nRows] = txt_segment(txt)
        %     H = histogram(txt.position(txt.incontact,2),100);
        % 
        %     rw.HVyrange = linspace(H.BinLimits(1),H.BinLimits(2),100);  % Range of y-values from histogram
        % 
        %     minRowDist = 10;    % Minimum distance between rows [mm]
        % 
        %     B = 1/3*ones(3,1);
        %     rw.HVfilt = filter(B,1,H.Values);   % Filtered histogram
        % 
        %     % Find maxima in histogram
        %     [~,iM] = findpeaks(rw.HVfilt,'MinPeakDistance',minRowDist-1);
        % 
        %     % Find minima in histogram
        %     [m,im] = findpeaks(-rw.HVfilt,'MinPeakDistance',minRowDist-1);
        % 
        %     rw.im = [1 im 100]; % Indx of minima in histogram
        %     rw.m = [0 m 0]; % Values of minima in histogram
        % 
        %     % Number of rows is denoted by number of maxima in histogram
        %     nRows = length(iM);
        % 
        %     for i=length(rw.im):-1:2
        %         % Find text points in every single row
        %         indxrow = find(txt.position(:,2)>rw.HVyrange(rw.im(i-1)) & txt.position(:,2)<rw.HVyrange(rw.im(i)));
        % 
        %         start_row = min(indxrow);   % Start of a single row
        %         end_row = max(indxrow); % End of a single row
        % 
        %         rows{i-1} = hw_row(txt,start_row,end_row);  % Call the hw_row class constructor for every single row
        %     end
        % end

        function [chunks,nChunks] = chunks_txt_segment(txt)
            figure
            plot_trajectory_incontact(txt)
            title('Chunks selection')

            n = 1;
            while 1
                point{n} = ginput;  % Select area vertices

                i = 1;
                while i<length(point{n})    % Draw selected area
                    hold on
                    plot([point{n}(i,1) point{n}(i+1,1)],[point{n}(i,2) point{n}(i+1,2)],'k')
                    i = i+1;
                end
                plot([point{n}(length(point{n}),1) point{n}(1,1)],[point{n}(length(point{n}),2) point{n}(1,2)],'k')

                % Find if points are inside the selected area
                isInside{n} = inpolygon(txt.position(:,1),txt.position(:,2),point{n}(:,1),point{n}(:,2));

                chunks{n} = hw_chunk(txt,isInside{n});  % Call the hw_chunk class constructor for every single chunks

                % To continue or not?
                answer = questdlg('Do you want to continue?','Choose an answer','Yes','No','No');

                switch answer
                    case 'Yes'
                        n = n+1;
                    case 'No'
                        break;
                end
            end

            nChunks = n;
        end

        function [words,nWords] = words_txt_segment(txt)
            figure(2)
            plot_trajectory_incontact(txt)
            title('Words selection')

            n = 1;
            while 1
                point{n} = ginput;
                
                i = 1;
                while i<length(point{n})
                    hold on
                    plot([point{n}(i,1) point{n}(i+1,1)],[point{n}(i,2) point{n}(i+1,2)],'k')
                    i = i+1;
                end
                figure(2)
                plot([point{n}(length(point{n}),1) point{n}(1,1)],[point{n}(length(point{n}),2) point{n}(1,2)],'k')

                isInside{n} = inpolygon(txt.position(:,1),txt.position(:,2),point{n}(:,1),point{n}(:,2));

                words{n} = hw_chunk(txt,isInside{n});   % Call the hw_chunk class constructor for every single word

                answer = questdlg('Do you want to continue?','Choose an answer','Yes','No','No');

                switch answer
                    case 'Yes'
                        n = n+1;
                    case 'No'
                        break;
                end
            end

            nWords = n;
        end

        function [theta,rho] = words_polar_coord(txt,wd)
            for i=1:txt.nWords
                % Words midpoints
                Xm(i) = mean(txt.words{i}.position(txt.words{i}.incontact,1));
                Ym(i) = mean(txt.words{i}.position(txt.words{i}.incontact,2));

                % Bring them to the reference point
                X(i) = Xm(i)-wd.Pm(1);
                Y(i) = Ym(i)-wd.Pm(2);

                % Polar coordinates
                [theta(i),rho(i)] = cart2pol(X(i),Y(i));
            end
        end

        function plot_words_polar_distribution(txt)
            polarscatter(txt.theta,txt.rho,'filled')
            title('Words distribution in text')
        end
    end
end