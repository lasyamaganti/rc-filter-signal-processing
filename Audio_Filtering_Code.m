%% ESE 3510 Case Study 1
% Authors: Ryan Soong r.l.soong@wustl.edu
%          Lasya Maganti l.maganti@wustl.edu

%frequency and angular frequency bands, manually tuned
close all;
clear
f1 = 500;
tau1 = 1/(2*pi*f1);

f2 = 2000; %low-mid
tau2 = 1/(2*pi*f2);

f3=3000; %some birds
tau3 = 1/(2*pi*f3);

f4 = 8000; %bird max
tau4 = 1/(2*pi*f4);

f5 = 15000; %very high
tau5 = 1/(2*pi*f5);

% Band1 -> Lowpass filter H_LP(s) = (1/tau) / (s + 1/tau)
b_band1 = [1/tau1];
a_band1= [1,1/tau1];

%Band2 -> Bandpass, lowpass at f2. Highpass at f1
% H_HP = s/(s+1/tau1)
%H_LP = (1/tau2)/(s+1/tau2)
b_hp_band2 = [1,0];
a_hp_band2 = [1,1/tau1];
b_lp_band2 = [1/tau2];
a_lp_band2 = [1,1/tau2];
%H = H_HP * H_LP
b_band2 = conv(b_hp_band2, b_lp_band2);
a_band2 = conv(a_hp_band2, a_lp_band2);

%Band3 -> bandpass, lowpass at f3, highpass at f2
b_hp_band3 = [1,0];
a_hp_band3 = [1,1/tau2];
b_lp_band3 = [1/tau3];
a_lp_band3 = [1,1/tau3];
%H = H_HP * H_LP
b_band3 = conv(b_hp_band3, b_lp_band3);
a_band3 = conv(a_hp_band3, a_lp_band3);

%Band4 -> bandpass, lowpass at f4, highpass at f3
b_hp_band4 = [1,0];
a_hp_band4 = [1,1/tau3];
b_lp_band4 = [1/tau4];
a_lp_band4 = [1,1/tau4];
%H = H_HP * H_LP
b_band4 = conv(b_hp_band4, b_lp_band4);
a_band4 = conv(a_hp_band4, a_lp_band4);

%Band 5 -> highpass at f5
% H_HP = s/(s+1/tau5)
b_band5= [1,0];
a_band5 = [1,1/tau5];

[xv,xvfs] = audioread('SNR Recording 2026-02-15 08_58.wav');
xv = xv(:,1);
t = [0:length(xv)-1]*1/xvfs;
%define gains
%magnify bird frequencies (band 4)
gain1= 1; gain2 = 1; gain3 = 1; gain4 = 1; gain5=1;

%process audio
y1 = lsim(b_band1, a_band1, xv, t);
y2 = lsim(b_band2, a_band2, xv, t);
y3 = lsim(b_band3, a_band3, xv, t);
y4 = lsim(b_band4, a_band4, xv, t);
y5 = lsim(b_band5, a_band5, xv, t);

y_new = (gain1*y1) + (gain2*y2) + (gain3*y3) + (gain4*y4) + (gain5*y5);
%normalize
y_new = y_new / max(abs(y_new));
sound(y_new, xvfs);


%% Visualization
%spectograms
figure;
hold on;
sgtitle('Spectrogram Plots: Unfiltered vs Filtered');
subplot(2,1,1);
xnormalised = xv/max(abs(xv));

spectrogram(xnormalised, 1024, 800, 1024, xvfs, 'yaxis'); %given
clim([-100 -80])
ylim([1 5])
title('Original Audio');
subplot(2,1,2);
spectrogram(y_new, 1024, 800, 1024, xvfs, 'yaxis')
clim([-80 -67])
ylim([1 5])
title('Filtered Audio');
hold off;

%% Bode plots, from hw3
frequencies = logspace(1,log10(20000),500);
w_bode = 2*pi*frequencies;
H1=freqs(b_band1,a_band1,w_bode); %talk about in report how this was much 
%more efficeint than complex exp. sweep because it wouldve required much
%more computation
H2=freqs(b_band2,a_band2,w_bode);
H3=freqs(b_band3,a_band3,w_bode);
H4=freqs(b_band4,a_band4,w_bode);
H5=freqs(b_band5,a_band5,w_bode);

H_bode_total = gain1*H1 + gain2*H2+ gain3*H3+ gain4*H4+ gain5*H5;
figure;

sgtitle('Equalizer Bode Plot');

%Magnitude plot
subplot(2,1,1);
semilogx(frequencies, 20*log10(abs(H1))); 
hold on;
semilogx(frequencies, 20*log10(abs(H2)));
semilogx(frequencies, 20*log10(abs(H3)));
semilogx(frequencies, 20*log10(abs(H4)));
semilogx(frequencies, 20*log10(abs(H5)));
semilogx(frequencies, 20*log10(abs(H_bode_total)), 'k--'); 
hold off;

title('Magnitude'); 
xlabel('Frequency (Hz)'); 
ylabel('Magnitude (dB)');

legend('B1', 'B2', 'B3', 'B4', 'B5', 'Total', 'Location', 'southwest');
grid on;
xlim([10 10^4*2])

%phase plot
subplot(2,1,2);
semilogx(frequencies, angle(H1)/pi); 
hold on;
semilogx(frequencies, angle(H2)/pi);
semilogx(frequencies, angle(H3)/pi);
semilogx(frequencies, angle(H4)/pi);
semilogx(frequencies, angle(H5)/pi);
semilogx(frequencies, angle(H_bode_total)/pi, 'k--');
hold off;
title('Phase Response'); 
xlabel('Frequency (Hz)'); 
ylabel('Phase (\times \pi rad)'); 
grid on;
xlim([10 10^4*2])


%impulse response graph
impulse = zeros(size(t));
impulse(1) = 1; %impulse at t=0

h1_impulse = lsim(b_band1, a_band1, impulse, t);
h2_impulse = lsim(b_band2, a_band2, impulse, t);
h3_impulse = lsim(b_band3, a_band3, impulse, t);
h4_impulse = lsim(b_band4, a_band4, impulse, t);
h5_impulse = lsim(b_band5, a_band5, impulse, t);

figure;
sgtitle('Impulse Response for Each Bands');
plot(t, h1_impulse);
hold on;
plot(t, h2_impulse);
plot(t, h3_impulse);
plot(t, h4_impulse);
plot(t, h5_impulse);
hold off;
xlabel('Time (s)');
ylabel('Amplitude');
legend('Band 1', 'Band 2', 'Band 3', 'Band 4', 'Band 5');
xlim([-0.00005 0.0003]); %to see begginging
grid on;

%time domain plot
figure;
subplot(2,1,1);
hold on;
plot(t,xv, 'Color','b');
title('Unfiltered');
subplot(2,1,2)

plot(t,y_new, 'Color', 'r');
title('Filtered');
hold off;
sgtitle('Time Domain: Unfiltered vs Filtered');
xlabel('Time (s)');
ylabel('Amplitude');



%Treble -> set gain4 and gain 5 to 4, others 0.3
% Bass -> gain1 and gain2 to 4, others 0.3
%unity -> all gains 1. 
%produce graphs for all these
%% Process audio file -> Space Station
[space_music, spacefs_music] = audioread('Space Station - Treble Cut.wav');
t_music = [0:length(space_music)-1] * 1/spacefs_music;
space_music = space_music(:,1);
% 5 bands
ym1 = lsim(b_band1, a_band1, space_music, t_music);
ym2 = lsim(b_band2, a_band2, space_music, t_music);
ym3 = lsim(b_band3, a_band3, space_music, t_music);
ym4 = lsim(b_band4, a_band4, space_music, t_music);
ym5 = lsim(b_band5, a_band5, space_music, t_music);

% Space Station recording
y_unity_space  = (1*ym1) + (1*ym2) + (1*ym3) + (1*ym4) + (1*ym5);
y_bass_space   = (4*ym1) + (4*ym2) + (0.3*ym3) + (0.3*ym4) + (0.3*ym5);
y_treble_space = (0.3*ym1) + (0.3*ym2) + (0.3*ym3) + (5*ym4) + (5*ym5);

y_unity_space  = y_unity_space  / max(abs(y_unity_space));
y_bass_space   = y_bass_space   / max(abs(y_bass_space));
y_treble_space = y_treble_space / max(abs(y_treble_space));
 
sound(y_unity_space, spacefs_music);

figure;
sgtitle('Time Domain: Space Station');

subplot(4,1,1);
plot(t_music, space_music, 'k');
title('Original Audio');

subplot(4,1,2);
plot(t_music, y_unity_space);
title('Unity Preset');

subplot(4,1,3);
plot(t_music, y_bass_space);
title('Bass Boost');

subplot(4,1,4);
plot(t_music, y_treble_space);
title('Treble Boost');
xlabel('Time (s)');

%spectograms
figure;
sgtitle('Spectograms: Space Station');

subplot(4,1,1);
spectrogram(space_music, 1024, 800, 1024, spacefs_music, 'yaxis');
ylabel(''); 
xlabel('');
clim([-80 -20])
ylim([0 5])
title('Original Audio');

subplot(4,1,2);
spectrogram(y_unity_space, 1024, 800, 1024, spacefs_music, 'yaxis');
ylabel(''); 
xlabel('');
title('Unity Preset');
clim([-80 -20])
ylim([0 5])


subplot(4,1,3);
spectrogram(y_bass_space, 1024, 800, 1024, spacefs_music, 'yaxis');
ylabel(''); 
xlabel('');
title('Bass Boost');
clim([-80 -20])
ylim([0 5])


subplot(4,1,4);
spectrogram(y_treble_space, 1024, 800, 1024, spacefs_music, 'yaxis');
ylabel(''); 
xlabel('');
title('Treble Boost');
clim([-80 -20])
ylim([0 5])


%% Process audio file -> Giant Steps
[giant_music, giantfs_music] = audioread('Giant Steps Bass Cut.wav');
t_music = [0:length(giant_music)-1] * 1/giantfs_music;
giant_music = giant_music(:,1);
% 5 bands
ym1 = lsim(b_band1, a_band1, giant_music, t_music);
ym2 = lsim(b_band2, a_band2, giant_music, t_music);
ym3 = lsim(b_band3, a_band3, giant_music, t_music);
ym4 = lsim(b_band4, a_band4, giant_music, t_music);
ym5 = lsim(b_band5, a_band5, giant_music, t_music);
%
% Giant Steps Recording
y_unity_giant  = (1*ym1) + (1*ym2) + (1*ym3) + (1*ym4) + (1*ym5);
y_bass_giant   = (4*ym1) + (4*ym2) + (0.3*ym3) + (0.3*ym4) + (0.3*ym5);
y_treble_giant = (0.3*ym1) + (0.3*ym2) + (0.3*ym3) + (5*ym4) + (5*ym5);

y_unity_giant  = y_unity_giant  / max(abs(y_unity_giant));
y_bass_giant   = y_bass_giant   / max(abs(y_bass_giant));
y_treble_giant = y_treble_giant / max(abs(y_treble_giant));

sound(y_unity_giant, giantfs_music);

figure;
sgtitle('Time Domain: Space Station');

subplot(4,1,1);
plot(t_music, giant_music, 'k');
title('Original Audio');

subplot(4,1,2);
plot(t_music, y_unity_giant);
title('Unity Preset');

subplot(4,1,3);
plot(t_music, y_bass_giant);
title('Bass Boost');

subplot(4,1,4);
plot(t_music, y_treble_giant);
title('Treble Boost');
xlabel('Time (s)');

%spectograms
figure;
sgtitle('Spectograms: Space Station');

subplot(4,1,1);
spectrogram(giant_music, 1024, 800, 1024, giantfs_music, 'yaxis');
ylabel(''); 
xlabel('');
clim([-80 -35])
ylim([0 5])
title('Original Audio');

subplot(4,1,2);
spectrogram(y_unity_giant, 1024, 800, 1024, giantfs_music, 'yaxis');
ylabel(''); 
xlabel('');
title('Unity Preset');
clim([-80 -35])
ylim([0 5])

subplot(4,1,3);
spectrogram(y_bass_giant, 1024, 800, 1024, giantfs_music, 'yaxis');
ylabel(''); 
xlabel('');
title('Bass Boost');
clim([-80 -35])
ylim([0 5])

subplot(4,1,4);
spectrogram(y_treble_giant, 1024, 800, 1024, giantfs_music, 'yaxis');
ylabel(''); 
xlabel('');
title('Treble Boost');
clim([-80 -35])
ylim([0 5])

audiowrite('bird_filtered.wav', y_new, xvfs);

audiowrite('space_unity.wav',  y_unity_space,  spacefs_music);
audiowrite('space_bass.wav',   y_bass_space,   spacefs_music);
audiowrite('space_treble.wav', y_treble_space, spacefs_music);

audiowrite('giant_unity.wav',  y_unity_giant,  giantfs_music);
audiowrite('giant_bass.wav',   y_bass_giant,   giantfs_music);
audiowrite('giant_treble.wav', y_treble_giant, giantfs_music);