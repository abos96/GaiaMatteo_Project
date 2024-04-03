classdef hw_row
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
        function rows = hw_row(txt,start_row,end_row)
            rows.time = txt.time(start_row:end_row);

            rows.position = txt.position(start_row:end_row,:);

            rows.pressure = txt.pressure(start_row:end_row);

            rows.incontact = find(rows.pressure>0);

            rows.inair = find(rows.pressure==0);

            [rows.pres_mean,rows.pres_std] = pressure_mean_std(rows);

            [rows.h,rows.ymin,rows.ymax] = height(rows);    % Calculate height

            [rows.l,rows.xmin,rows.xmax] = length(rows);    % Calculate length

            [rows.vel,rows.speed,rows.speed_mean,rows.speed_std] = velocity(rows,txt);

            [rows.acc,rows.acc_mean,rows.acc_std] = acceleration(rows,txt);

            [rows.jer,rows.jer_norm,rows.jer_mean,rows.jer_std] = jerk(rows,txt);

            [rows.Rc,rows.Rc_mean,rows.Rc_std] = curvrad(rows);

            [rows.Pc,rows.Pa] = time_percentage(rows);

            rows.t = tilt(rows);    % Calculate tilt
        end

        function [pres_mean,pres_std] = pressure_mean_std(rows)
            pres_mean = mean(rows.pressure(rows.incontact));
            pres_std = std(rows.pressure(rows.incontact));
        end

        function plot_pressure(rows,label)
            line(rows.time,rows.pressure,'col','b')
            line([rows.time(1) rows.time(end)],[rows.pres_mean rows.pres_mean],'col','r')
            xlabel('Time [s]')
            ylabel('P [MPa]')
            title(label+' pressure')
            legend('Pressure','Mean pressure')
        end

        function plot_pressure_distribution(rows,label)
            pressure_no_outliers = rmoutliers(rows.pressure);

            boxplot(pressure_no_outliers)
            title(label+' pressure distribution plot')
            xticklabels('Pressure')
        end

        function [h,ymin,ymax] = height(rows)
            ymax = median(findpeaks(rows.position(rows.incontact,2)));  % Median of maximum peaks on y
            ymin = median(-findpeaks(-rows.position(rows.incontact,2)));    % Median of minimum peaks on y
            h = ymax-ymin;
        end

        function [l,xmin,xmax] = length(rows)
            xmax = max(rows.position(rows.incontact,1));    % Maximum position on x
            xmin = min(rows.position(rows.incontact,1));    % Minimum position on x
            l = xmax-xmin;
        end

        function [vel,speed,speed_mean,speed_std] = velocity(rows,txt)
            vx = ((-txt.fs)^1)*savitzkyGolayFilt(rows.position(:,1),txt.SGorder,1,txt.SGwin);
            vy = ((-txt.fs)^1)*savitzkyGolayFilt(rows.position(:,2),txt.SGorder,1,txt.SGwin);

            vel.incontact = nan(length(rows.time),2);
            vel.incontact(rows.incontact,:) = [vx(rows.incontact) vy(rows.incontact)];

            vel.inair = nan(length(rows.time),2);
            vel.inair(rows.inair,:) = [vx(rows.inair) vy(rows.inair)];

            speed.incontact = nan(length(rows.time),1);
            speed.incontact(rows.incontact) = sqrt(vx(rows.incontact).^2+vy(rows.incontact).^2);

            speed.inair = nan(length(rows.time),1);
            speed.inair(rows.inair) = sqrt(vx(rows.inair).^2+vy(rows.inair).^2);

            speed_mean.incontact = nanmean(speed.incontact);
            speed_std.incontact = nanstd(speed.incontact);

            speed_mean.inair = nanmean(speed.inair);
            speed_std.inair = nanstd(speed.inair);
        end

        function plot_speed_incontact(rows,label)
            line(rows.time,rows.speed.incontact)
            xlabel('Time [s]')
            ylabel('V [mm/s]')
            title(label+' speed in contact')
        end

        function plot_speed_inair(rows,label)
            line(rows.time,rows.speed.inair)
            xlabel('Time [s]')
            ylabel('V [mm/s]')
            title(label+' speed in air')
        end

        function plot_speed_distribution(rows,label)
            speed_no_outliers = rmoutliers([rows.speed.incontact rows.speed.inair]);

            boxplot(speed_no_outliers)
            title(label+' speed distribution plot')
            xticklabels({'Speed in contact','Speed in air'})
        end

        function [acc,acc_mean,acc_std] = acceleration(rows,txt)
            ax = ((-txt.fs)^2)*savitzkyGolayFilt(rows.position(:,1),txt.SGorder,2,txt.SGwin);
            ay = ((-txt.fs)^2)*savitzkyGolayFilt(rows.position(:,2),txt.SGorder,2,txt.SGwin);

            acc.incontact = nan(length(rows.time),2);
            acc.incontact(rows.incontact,:) = [ax(rows.incontact) ay(rows.incontact)];

            acc.inair = nan(length(rows.time),2);
            acc.inair(rows.inair,:) = [ax(rows.inair) ay(rows.inair)];

            acc_mean.incontact = nanmean(acc.incontact);
            acc_std.incontact = nanstd(acc.incontact);

            acc_mean.inair = nanmean(acc.inair);
            acc_std.inair = nanstd(acc.inair);
        end

        function plot_acceleration_incontact(rows,label)
            subplot(2,1,1)
            line(rows.time,rows.acc.incontact(:,1),'col','b')
            ylabel('A_x [mm/s^2]')
            title(label+' acceleration in contact')
            subplot(2,1,2)
            line(rows.time,rows.acc.incontact(:,2),'col','b')
            xlabel('Time [s]')
            ylabel('A_y [mm/s^2]')
        end

        function plot_acceleration_inair(rows,label)
            subplot(2,1,1)
            line(rows.time,rows.acc.inair(:,1),'col','b')
            ylabel('A_x [mm/s^2]')
            title(label+' acceleration in air')
            subplot(2,1,2)
            line(rows.time,rows.acc.inair(:,2),'col','b')
            xlabel('Time [s]')
            ylabel('A_y [mm/s^2]')
        end

        function plot_acceleration_distribution(rows,label)
            acc_no_outliers = rmoutliers([rows.acc.incontact rows.acc.inair]);

            boxplot(acc_no_outliers)
            title(label+' acceleration distribution plot')
            xticklabels({'Ax in contact','Ay in contact','Ax in air','Ay in air'})
        end

        function [jer,jer_norm,jer_mean,jer_std] = jerk(rows,txt)
            jx = ((-txt.fs)^3)*savitzkyGolayFilt(rows.position(:,1),txt.SGorder,3,txt.SGwin);
            jy = ((-txt.fs)^3)*savitzkyGolayFilt(rows.position(:,2),txt.SGorder,3,txt.SGwin);

            jer.incontact = nan(length(rows.time),2);
            jer.incontact(rows.incontact,:) = [jx(rows.incontact) jy(rows.incontact)];

            jer.inair = nan(length(rows.time),2);
            jer.inair(rows.inair,:) = [jx(rows.inair) jy(rows.inair)];

            jer_norm.incontact = nan(length(rows.time),1);
            jer_norm.incontact(rows.incontact) = sqrt(jx(rows.incontact).^2+jy(rows.incontact).^2);

            jer_norm.inair = nan(length(rows.time),1);
            jer_norm.inair(rows.inair) = sqrt(jx(rows.inair).^2+jy(rows.inair).^2);

            jer_mean.incontact = nanmean(jer_norm.incontact);
            jer_std.incontact = nanstd(jer_norm.incontact);

            jer_mean.inair = nanmean(jer_norm.inair);
            jer_std.inair = nanstd(jer_norm.inair);
        end

        function plot_jerk_incontact(rows,label)
            line(rows.time,rows.jer_norm.incontact)
            xlabel('Time [s]')
            ylabel('J [mm/s^3]')
            title(label+' jerk in contact')
        end

        function plot_jerk_inair(rows,label)
            line(rows.time,rows.jer_norm.inair)
            xlabel('Time [s]')
            ylabel('J [mm/s^3]')
            title(label+' jerk in air')
        end

        function plot_jerk_distribution(rows,label)
            jerk_no_outliers = rmoutliers([rows.jer_norm.incontact rows.jer_norm.inair]);

            boxplot(jerk_no_outliers)
            title(label+' jerk distribution plot')
            xticklabels({'Jerk in contact','Jerk in air'})
        end

        function [Rc,Rc_mean,Rc_std] = curvrad(rows)
            k = abs(rows.vel.incontact(:,1).*rows.acc.incontact(:,2)-rows.vel.incontact(:,2).*rows.acc.incontact(:,1))./((rows.vel.incontact(:,1).^2+rows.vel.incontact(:,2).^2).^(3/2));
            Rc = 1./k;

            Rc_mean = nanmean(Rc);
            Rc_std = nanstd(Rc);
        end

        function plot_curvrad(rows,label)
            line(rows.time,rows.Rc)
            xlabel('Time [s]')
            ylabel('Rc [mm]')
            title(label+' curvature radius')
            set(gca,'YScale','log')
        end

        function plot_curvrad_distribution(rows,label)
            curvrad_no_outliers = rmoutliers(rows.Rc);

            boxplot(curvrad_no_outliers)
            title(label+' curvature radius distribution plot')
            xticklabels('Curvature radius')
        end

        function [Pc,Pa] = time_percentage(rows)
            Pc = length(rows.time(rows.incontact))/length(rows.time)*100;
            Pa = length(rows.time(rows.inair))/length(rows.time)*100;
        end

        function t = tilt(rows)
            % Tilt calculation from pca
            U = pca(rows.position(rows.incontact,:));

            if U(1,1)<0
                U(1,1) = -U(1,1);
            end

            t = atan2d(U(1,2),U(1,1));  % Calculate atan2 (four-quadrant inverse tangent) in degrees
        end
    end
end