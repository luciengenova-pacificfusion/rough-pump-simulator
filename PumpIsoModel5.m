% Variables (all units SI)
clear; clc; close all;
pkg load control
pkg load signal

% ----- Isolator selection -----
% Set to 1, 2, or 3 to choose which isolator is processed and plotted:
% 1 = Sorbothane, 2 = Air shock, 3 = Wire rope
% Set to 4 to process and plot all three isolators
isoSelect = 4;
% ------------------------------

Mpump = 75;
Mheatex = 7;
Mbase = 22*2;

Mdpm = 16000;
Mairbox = 4000;

Mopt = 877.8/2.2 % Mass of optical table top surface
KoptLeg = 1766396.89; % Ns/m
Kopt = KoptLeg*4 % Spring rate of all 4 optical table legs
ZoptLeg = 0.18265;
CoptLeg = 2 * ZoptLeg * sqrt(KoptLeg * Mopt);
Copt = CoptLeg*4 % Damping coefficient of all 4 optical table legs


% ------ Setup Selection ----
% Set to 1 to simulate lab optical table test setup
% Set to 2 to simulate DPM airbox only
% Set to 3 or other to simulate DPM
setup = 2;

m1 = Mpump + Mbase + Mheatex % kg
if setup == 1
  m2 = Mopt;
elseif setup == 2
  m2 = Mairbox;
else
  m2 = Mdpm;
end

F1 = 60; % Estimated force amplitude (N)
f1 = 112; % forcing frequency (Hz)
w1 = 2*pi*f1; % forcing frequency (rad/s)

F2 = 20; % Estimated force amplitude (N)
f2 = 62; % forcing frequency (Hz)
w2 = 2*pi*f2; % forcing frequency (rad/s)

F3 = 60; % Estimated  force amplitude (N)
f3 = 224; % forcing frequency (Hz)
w3 = 2*pi*f3; % forcing frequency (rad/s)

isolators = {'Sorbothane','Air shock','Wire rope'};

% Resolve selection into the set of isolator indices to process/plot
if isoSelect == 4
    isoRange = 1:length(isolators);
else
    isoRange = isoSelect;
end

##k1_list = [166706.80*6, 35110.95*4, 101281.59*4];   % Researched Values N/m
##k1_list = [83829*6, 33125*4, 472333*4];   % Static Deflection Values N/m
k1_list = [9256*6, 33125*4, 47233*4];   % Match Tuned Values N/m
##z1_list = [0.18265*6, 0.05*4, 0.2*4]; % Researched Values Ns/m
z1_list = [0.18*6, 0.1*4, 0.2*4]; % Match Tuned Values Ns/m

wn_list = sqrt(k1_list./m1)/(2*pi) % Natural frequencies in Hz

% Estimate damping coefficients based on damping ratio, mass, and stiffness
% Single mass model assumption valid as m1 << m2 (meff = m1m2/(m1+m2))
c1_list = 2 .* z1_list .* sqrt(k1_list .* m1);    % Ns/m
k2 = Kopt;
c2 = Copt;
##k2 = 0; % Worst case model
##c2 = 0; % Worst case model

% Vibration specs
VrmsToXmag = 1.414*(1/w1);
vccMag = (12.5e-6)*VrmsToXmag;
vcdMag = (6.25e-6)*VrmsToXmag;

% VC limits
vccVrms = 12.5e-6;          % m/s RMS
vccVelocityPeak = sqrt(2)*vccVrms;  % m/s peak form

% Time response
s = tf('s');

for i = isoRange

    k1 = k1_list(i);
    c1 = c1_list(i);

    % Response to forcing
    len = 10; % Sample in seconds
    resolution = 10000; % Samples per second
    samples = len*resolution;

    t = linspace(0,len,samples);
    A = F1*sin(w1*t) + F2*sin(w2*t) + F3*sin(w3*t);

    % Set up transfer function (output/input)
    x2x1 = (c1*s + k1)/(m2*s^2 + (c1+c2)*s + (k1+k2));
    x2A = (c1*s + k1)/((m1*s^2 + c1*s + k1)*(m2*s^2 + (c1+c2)*s + (k1+k2)) - (c1*s + k1)^2); % Pump forcing input, DPM displacement output
    x2A_array{i} = x2A;
    x1A = x2A / x2x1;
    x1A_array{i} = x1A;

    a2A = s^2*x2A;      % acceleration transfer function
    a2A_array{i} = a2A;
    a1A = a2A / x2x1;
    a1A_array{i} = a1A;

    v2A = s*x2A;      % velocity transfer function
    v2A_array{i} = v2A;
    v1A = v2A / x2x1;
    v1A_array{i} = v1A;

    % Gain at operating frequency
    [mag, phase] = bode(x2A, w1);
    gain_dB = mag2db(mag);
##    pIso = (1 - mag(:)/1)*100

end

% Time response acceleration
figure(2);

colors = lines(length(isolators));

% Pump acceleration
subplot(2,1,1);
hold on
grid on
for i = isoRange
    acc1 = lsim(a1A_array{i}, A, t); % pump accel
    plot(t, acc1, 'Color', colors(i,:), 'DisplayName', isolators{i});
end
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Pump Acceleration Time Response')
legend('Location','northeast')
hold off

% Base acceleration
subplot(2,1,2);
hold on
grid on
for i = isoRange
    acc2 = lsim(a2A_array{i}, A, t); % base accel
    plot(t, acc2, 'Color', colors(i,:), 'DisplayName', isolators{i});
end
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Base Acceleration Time Response')
legend('Location','northeast')
hold off

% Welch PSD - RMS velocity spectrum
figure(5);

fs = resolution;                 % sampling frequency (Hz)
nfft = 2^15;                     % FFT length
window = hanning(nfft);
noverlap = 0.5;

% Enable coherent scaling of Welch PSD to RMS per bin
% RMS velocity per bin = sqrt(PSD * df)
df = fs/nfft;

% VC limit lines (RMS velocity, um/s)
vcC = 12.5;   % VC-C
vcD = 6.25;   % VC-D
vcE = 3.12;   % VC-E

% Pump velocity spectrum
subplot(2,1,1);
hold on
grid on
for i = isoRange
    vel1 = lsim(v1A_array{i}, A, t);   % pump velocity time response
    [Pxx, fw] = pwelch(vel1, window, noverlap, nfft, fs);
    vrms = sqrt(Pxx * df)*1e6;
    plot(fw, vrms, 'DisplayName', isolators{i});
end
xlabel('Frequency (Hz)')
ylabel('RMS Velocity (um/s)')
title('Pump Welch RMS Velocity Spectrum')
set(gca, 'XScale', 'linear', 'YScale', 'linear');
xlim([0 1000]);
xticks(0:50:1000);
yline(vcC, '--', 'VC-C');
yline(vcD, '--', 'VC-D');
yline(vcE, '--', 'VC-E');
##legend('Location','northeast')
hold off

% Base velocity spectrum
subplot(2,1,2);
hold on
grid on
for i = isoRange
    vel2 = lsim(v2A_array{i}, A, t);   % base velocity time response
    [Pxx, fw] = pwelch(vel2, window, noverlap, nfft, fs);
    vrms = sqrt(Pxx * df)*1e6;
    plot(fw, vrms, 'DisplayName', isolators{i});
end
xlabel('Frequency (Hz)')
ylabel('RMS Velocity (um/s)')
title('Base Welch RMS Velocity Spectrum')
set(gca, 'XScale', 'linear', 'YScale', 'linear');
xlim([0 1000]);
xticks(0:50:1000);
yline(vcC, '--', 'VC-C');
yline(vcD, '--', 'VC-D');
yline(vcE, '--', 'VC-E');
##legend('Location','northeast')
hold off

##% Time response position
##figure(1);
##hold on
##grid on
##for i = isoRange
##    iso = lsim(x2A_array{i}, A, t); % base accel
##    plot(t, iso, 'Color', colors(i,:), 'DisplayName', isolators{i});
##end
####fprintf('\n%s\n', isolators{i});
####fprintf('k1 = %.2f N/m, c1 = %.2f Ns/m\n', k1, c1);
####fprintf('Operating Freq = %.3f Hz (%.2f rad/s)\n', f, w);
####fprintf('Gain = %.3e (%.2f dB), Phase = %.2f deg\n', mag(:), gain_dB(:), phase(:));
####fprintf('\n');
##xlabel('Time (s)')
##ylabel('Displacement (m)')
##title('Displacement Time Response')
##legend('Location','northeast')
####% VC limit markers
####yline(vccMag, '--', 'VC-C', 'DisplayName', 'VC-C');
####yline(-vccMag, '--', 'HandleVisibility','off');
##hold off

##% Time response velocity
##figure(3);
##hold on
##grid on
##
##for i = isoRange
####    acc = lsim(v1A_array{i}, A, t); % plot pump accel
##    acc = lsim(v2A_array{i}, A, t); % plot base accel
##    plot(t, acc, 'DisplayName', isolators{i});
##end
##
##xlabel('Time (s)')
##ylabel('Velocity (m/s)')
##title('Velocity Time Response')
##legend('Location','northeast')
##
##% VC limit markers
##yline(vccVelocityPeak,'--','VC-C','DisplayName','VC-C');
##yline(-vccVelocityPeak,'--','HandleVisibility','off');
##
##hold off
##%%

##% Frequency response
##figure(4); clf
##hold on
##
##w_plot = logspace(0,3,1000); % rad/s
##
##for i = isoRange
##    [mag, ~] = bode(x2A_array{i}, w_plot);
##    mag = squeeze(mag);
##    semilogx(w_plot/(2*pi), mag2db(mag), 'DisplayName', isolators{i});
##end
##
##grid on
##xlabel('Frequency (Hz)')
##ylabel('Magnitude (dB)')
##title('Frequency Response')
##legend('Location','northeast')
##
##% Operating frequency marker
##xline(f, '--k', sprintf('Operating Freq = %.2f Hz', f), 'LabelVerticalAlignment', 'middle', 'LabelOrientation', 'horizontal');
##
##hold off
