clear all;
close all;
clc;
%%
[file,pathdir] = uigetfile(fullfile('Prove','*.csv'));  % Select data file
orient = 'l';   % l o p

wd = read(wacomdata,fullfile(pathdir,file),orient); % Call the wacomdata class constructor
%%
% Savitzky-Golay parameters
SGorder = 4;
SGwin = 11;

txt = hw_text(wd,SGwin,SGorder);    % Call the hm_text class constructor
%%
% Distribution of words (polar coordinates)
figure
plot_words_polar_distribution(txt)
%%
figure
plot_trajectory(txt)
%%
% Show selected chunks
for i=1:txt.nChunks
    figure
    plot(txt.chunks{i}.position(txt.chunks{i}.incontact,1),txt.chunks{i}.position(txt.chunks{i}.incontact,2),'marker','.','lines','none')
    title("Chunk "+i)
    axis equal
end

% Show selected words
for i=1:txt.nWords
    figure
    plot(txt.words{i}.position(txt.words{i}.incontact,1),txt.words{i}.position(txt.words{i}.incontact,2),'marker','.','lines','none')
    title("Word "+i)
    axis equal
end
%% Text kinetic analysis
figure
plot_pressure(txt)

figure
plot_speed_incontact(txt)

figure
plot_speed_inair(txt)

figure
plot_acceleration_incontact(txt)

figure
plot_acceleration_inair(txt)

figure
plot_jerk_incontact(txt)

figure
plot_jerk_inair(txt)

figure
plot_curvrad(txt)
%% Text statistical analysis
figure
plot_pressure_distribution(txt)

figure
plot_speed_distribution(txt)

figure
plot_acceleration_distribution(txt)

figure
plot_jerk_distribution(txt)

figure
plot_curvrad_distribution(txt)
%% Rows kinetic analysis
% rows = txt.rows;    % Variable rename
% 
% [r,label] = select_row(txt);    % Select row number
% 
% figure
% plot_pressure(rows{r},label)
% 
% figure
% plot_speed_incontact(rows{r},label)
% 
% figure
% plot_speed_inair(rows{r},label)
% 
% figure
% plot_acceleration_incontact(rows{r},label)
% 
% figure
% plot_acceleration_inair(rows{r},label)
% 
% figure
% plot_jerk_incontact(rows{r},label)
% 
% figure
% plot_jerk_inair(rows{r},label)
% 
% figure
% plot_curvrad(rows{r},label)
%% Rows statistical analysis
% figure
% plot_pressure_distribution(rows{r},label)
% 
% figure
% plot_speed_distribution(rows{r},label)
% 
% figure
% plot_acceleration_distribution(rows{r},label)
% 
% figure
% plot_jerk_distribution(rows{r},label)
% 
% figure
% plot_curvrad_distribution(rows{r},label)
%% Chunks kinetic analysis
for i=1:txt.nChunks
    figure
    plot_pressure(txt.chunks{i},i)

    figure
    plot_speed_incontact(txt.chunks{i},i)

    figure
    plot_speed_inair(txt.chunks{i},i)

    figure
    plot_acceleration_incontact(txt.chunks{i},i)

    figure
    plot_acceleration_inair(txt.chunks{i},i)

    figure
    plot_jerk_incontact(txt.chunks{i},i)

    figure
    plot_jerk_inair(txt.chunks{i},i)

    figure
    plot_curvrad(txt.chunks{i},i)
end
%% Chunks statistical analysis
for i=1:txt.nChunks
    figure
    plot_pressure_distribution(txt.chunks{i},i)

    figure
    plot_speed_distribution(txt.chunks{i},i)

    figure
    plot_acceleration_distribution(txt.chunks{i},i)

    figure
    plot_jerk_distribution(txt.chunks{i},i)

    figure
    plot_curvrad_distribution(txt.chunks{i},i)
end
%% Words kinetic analysis
for i=1:txt.nWords
    figure
    plot_pressure(txt.words{i},i)

    figure
    plot_speed_incontact(txt.words{i},i)

    figure
    plot_speed_inair(txt.words{i},i)

    figure
    plot_acceleration_incontact(txt.words{i},i)

    figure
    plot_acceleration_inair(txt.words{i},i)

    figure
    plot_jerk_incontact(txt.words{i},i)

    figure
    plot_jerk_inair(txt.words{i},i)

    figure
    plot_curvrad(txt.words{i},i)
end
%% Words statistical analysis
for i=1:txt.nWords
    figure
    plot_pressure_distribution(txt.words{i},i)

    figure
    plot_speed_distribution(txt.words{i},i)

    figure
    plot_acceleration_distribution(txt.words{i},i)

    figure
    plot_jerk_distribution(txt.words{i},i)

    figure
    plot_curvrad_distribution(txt.words{i},i)
end