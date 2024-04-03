classdef hw_chunk
    properties
        time = [];
        position = [];
        pressure = [];
        incontact = [];
        inair = [];
        pres_mean = [];
        pres_std = [];
        h = [];
        ymin = [];
        ymax = [];
        l = [];
        xmin = [];
        xmax = [];
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
        t = [];
    end

    methods
        function chunk = hw_chunk(txt,isInside)
            chunk.time = txt.time(find(isInside));

            chunk.position = txt.position(find(isInside),:);

            chunk.pressure = txt.pressure(find(isInside));

            chunk.incontact = find(chunk.pressure>0);

            chunk.inair = find(chunk.pressure==0);

            [chunk.pres_mean,chunk.pres_std] = pressure_mean_std(chunk);

            [chunk.h,chunk.ymin,chunk.ymax] = height_chunk(chunk);    % Calculate height

            [chunk.l,chunk.xmin,chunk.xmax] = length_chunk(chunk);    % Calculate length

            [chunk.vel,chunk.speed,chunk.speed_mean,chunk.speed_std] = velocity(chunk,txt);

            [chunk.acc,chunk.acc_mean,chunk.acc_std] = acceleration(chunk,txt);

            [chunk.jer,chunk.jer_norm,chunk.jer_mean,chunk.jer_std] = jerk(chunk,txt);

            [chunk.Rc,chunk.Rc_mean,chunk.Rc_std] = curvrad(chunk);

            [chunk.Pc,chunk.Pa] = time_percentage(chunk);

            chunk.t = tilt(chunk);  % Calculate tilt
        end

        function [pres_mean,pres_std] = pressure_mean_std(chunk)
            pres_mean = mean(chunk.pressure(chunk.incontact));
            pres_std = std(chunk.pressure(chunk.incontact));
        end

        function plot_pressure(chunk,n)
            line(chunk.time,chunk.pressure,'col','b')
            line([chunk.time(1) chunk.time(end)],[chunk.pres_mean chunk.pres_mean],'col','r')
            xlabel('Time [s]')
            ylabel('P levels')
            title("Row "+n+" pressure")
            legend('Pressure','Mean pressure')
        end

        function plot_pressure_distribution(chunk,n)
            pressure_no_outliers = rmoutliers(chunk.pressure);

            boxplot(pressure_no_outliers)
            title("Row "+n+" pressure distribution plot")
            xticklabels('Pressure')
        end

        function [h,ymin,ymax] = height_chunk(chunk)
            ymax = mean(findpeaks(chunk.position(chunk.incontact,2)));  % Mean of maximum peaks on y
            ymin = mean(-findpeaks(-chunk.position(chunk.incontact,2)));    % Mean of minimum peaks on y
            h = ymax-ymin;
        end

        function [l,xmin,xmax] = length_chunk(chunk)
            xmax = max(chunk.position(chunk.incontact,1));  % Maximum position on x
            xmin = min(chunk.position(chunk.incontact,1));  % Minimum position on x
            l = xmax-xmin;
        end

        function [vel,speed,speed_mean,speed_std] = velocity(chunk,txt)
            vx = ((-txt.fs)^1)*savitzkyGolayFilt(chunk.position(:,1),txt.SGorder,1,txt.SGwin);
            vy = ((-txt.fs)^1)*savitzkyGolayFilt(chunk.position(:,2),txt.SGorder,1,txt.SGwin);

            vel.incontact = nan(length(chunk.time),2);
            vel.incontact(chunk.incontact,:) = [vx(chunk.incontact) vy(chunk.incontact)];

            vel.inair = nan(length(chunk.time),2);
            vel.inair(chunk.inair,:) = [vx(chunk.inair) vy(chunk.inair)];

            speed.incontact = nan(length(chunk.time),1);
            speed.incontact(chunk.incontact) = sqrt(vx(chunk.incontact).^2+vy(chunk.incontact).^2);

            speed.inair = nan(length(chunk.time),1);
            speed.inair(chunk.inair) = sqrt(vx(chunk.inair).^2+vy(chunk.inair).^2);

            speed_mean.incontact = mean(speed.incontact,"omitnan");
            speed_std.incontact = std(speed.incontact,"omitnan");

            speed_mean.inair = mean(speed.inair,"omitnan");
            speed_std.inair = std(speed.inair,"omitnan");
        end

        function plot_speed_incontact(chunk,n)
            line(chunk.time,chunk.speed.incontact)
            xlabel('Time [s]')
            ylabel('V [mm/s]')
            title("Row "+n+" speed in contact")
        end

        function plot_speed_inair(chunk,n)
            line(chunk.time,chunk.speed.inair)
            xlabel('Time [s]')
            ylabel('V [mm/s]')
            title("Row "+n+" speed in air")
        end

        function plot_speed_distribution(chunk,n)
            speed_no_outliers = rmoutliers([chunk.speed.incontact chunk.speed.inair]);

            boxplot(speed_no_outliers)
            title("Row "+n+" speed distribution plot")
            xticklabels({'Speed in contact','Speed in air'})
        end

        function [acc,acc_mean,acc_std] = acceleration(chunk,txt)
            ax = ((-txt.fs)^2)*savitzkyGolayFilt(chunk.position(:,1),txt.SGorder,2,txt.SGwin);
            ay = ((-txt.fs)^2)*savitzkyGolayFilt(chunk.position(:,2),txt.SGorder,2,txt.SGwin);

            acc.incontact = nan(length(chunk.time),2);
            acc.incontact(chunk.incontact,:) = [ax(chunk.incontact) ay(chunk.incontact)];

            acc.inair = nan(length(chunk.time),2);
            acc.inair(chunk.inair,:) = [ax(chunk.inair) ay(chunk.inair)];

            acc_mean.incontact = mean(acc.incontact,"omitnan");
            acc_std.incontact = std(acc.incontact,"omitnan");

            acc_mean.inair = mean(acc.inair,"omitnan");
            acc_std.inair = std(acc.inair,"omitnan");
        end

        function plot_acceleration_incontact(chunk,n)
            subplot(2,1,1)
            line(chunk.time,chunk.acc.incontact(:,1),'col','b')
            ylabel('A_x [mm/s^2]')
            title("Row "+n+" acceleration in contact")
            subplot(2,1,2)
            line(chunk.time,chunk.acc.incontact(:,2),'col','b')
            xlabel('Time [s]')
            ylabel('A_y [mm/s^2]')
        end

        function plot_acceleration_inair(chunk,n)
            subplot(2,1,1)
            line(chunk.time,chunk.acc.inair(:,1),'col','b')
            ylabel('A_x [mm/s^2]')
            title("Row "+n+" acceleration in air")
            subplot(2,1,2)
            line(chunk.time,chunk.acc.inair(:,2),'col','b')
            xlabel('Time [s]')
            ylabel('A_y [mm/s^2]')
        end

        function plot_acceleration_distribution(chunk,n)
            acc_no_outliers = rmoutliers([chunk.acc.incontact chunk.acc.inair]);

            boxplot(acc_no_outliers)
            title("Row "+n+" acceleration distribution plot")
            xticklabels({'Ax in contact','Ay in contact','Ax in air','Ay in air'})
        end

        function [jer,jer_norm,jer_mean,jer_std] = jerk(chunk,txt)
            jx = ((-txt.fs)^3)*savitzkyGolayFilt(chunk.position(:,1),txt.SGorder,3,txt.SGwin);
            jy = ((-txt.fs)^3)*savitzkyGolayFilt(chunk.position(:,2),txt.SGorder,3,txt.SGwin);

            jer.incontact = nan(length(chunk.time),2);
            jer.incontact(chunk.incontact,:) = [jx(chunk.incontact) jy(chunk.incontact)];

            jer.inair = nan(length(chunk.time),2);
            jer.inair(chunk.inair,:) = [jx(chunk.inair) jy(chunk.inair)];

            jer_norm.incontact = nan(length(chunk.time),1);
            jer_norm.incontact(chunk.incontact) = sqrt(jx(chunk.incontact).^2+jy(chunk.incontact).^2);

            jer_norm.inair = nan(length(chunk.time),1);
            jer_norm.inair(chunk.inair) = sqrt(jx(chunk.inair).^2+jy(chunk.inair).^2);

            jer_mean.incontact = mean(jer_norm.incontact,"omitnan");
            jer_std.incontact = std(jer_norm.incontact,"omitnan");

            jer_mean.inair = mean(jer_norm.inair,"omitnan");
            jer_std.inair = std(jer_norm.inair,"omitnan");
        end

        function plot_jerk_incontact(chunk,n)
            line(chunk.time,chunk.jer_norm.incontact)
            xlabel('Time [s]')
            ylabel('J [mm/s^3]')
            title("Row "+n+" jerk in contact")
        end

        function plot_jerk_inair(chunk,n)
            line(chunk.time,chunk.jer_norm.inair)
            xlabel('Time [s]')
            ylabel('J [mm/s^3]')
            title("Row "+n+" jerk in air")
        end

        function plot_jerk_distribution(chunk,n)
            jerk_no_outliers = rmoutliers([chunk.jer_norm.incontact chunk.jer_norm.inair]);

            boxplot(jerk_no_outliers)
            title("Row "+n+" jerk distribution plot")
            xticklabels({'Jerk in contact','Jerk in air'})
        end

        function [Rc,Rc_mean,Rc_std] = curvrad(chunk)
            k = abs(chunk.vel.incontact(:,1).*chunk.acc.incontact(:,2)-chunk.vel.incontact(:,2).*chunk.acc.incontact(:,1))./((chunk.vel.incontact(:,1).^2+chunk.vel.incontact(:,2).^2).^(3/2));
            Rc = 1./k;

            Rc_mean = mean(Rc,"omitnan");
            Rc_std = std(Rc,"omitnan");
        end

        function plot_curvrad(chunk,n)
            line(chunk.time,chunk.Rc)
            xlabel('Time [s]')
            ylabel('Rc [mm]')
            title("Row "+n+" curvature radius")
            set(gca,'YScale','log')
        end

        function plot_curvrad_distribution(chunk,n)
            curvrad_no_outliers = rmoutliers(chunk.Rc);

            boxplot(curvrad_no_outliers)
            title("Row "+n+" curvature radius distribution plot")
            xticklabels('Curvature radius')
        end

        function [Pc,Pa] = time_percentage(chunk)
            Pc = length(chunk.time(chunk.incontact))/length(chunk.time)*100;
            Pa = length(chunk.time(chunk.inair))/length(chunk.time)*100;
        end

        function t = tilt(chunk)
            % Tilt calculation from pca
            U = pca(chunk.position(chunk.incontact,:));

            if U(1,1)<0
                U(1,1) = -U(1,1);
            end

            t = atan2d(U(1,2),U(1,1));  % Calculate atan2 (four-quadrant inverse tangent) in degrees
        end
    end
end