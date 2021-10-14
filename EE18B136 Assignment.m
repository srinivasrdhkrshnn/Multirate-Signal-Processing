clc ;
clear all ;
close all ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
wp = 0.45*pi ;
ws = 0.55*pi ;
del = 1 ;
attenuation = 60 ;
ds_factor = 2 ;
us_factor = 3 ;
intfilt_support = 3 ;
alpha = 1 ;
interpolator = intfilt(us_factor,intfilt_support,alpha) ;
fvtool(interpolator) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Original Signal

[y_08,fs_08] = audioread('speech8khz.wav') ;
[y_16,fs_16] = audioread('music16khz.wav') ;

N_08 = length(y_08) ;
N_16 = length(y_16) ;

filter_08 = equiripple_filter(fs_08,wp,ws,del,attenuation);
filter_16 = equiripple_filter(fs_16,wp,ws,del,attenuation);
fvtool(filter_08) ;
% fvtool(filter_16) ;

FT_08 = abs(fftshift(fft(y_08,N_08))) ;
FT_16 = abs(fftshift(fft(y_16,N_16))) ;

f_08 = (0:N_08-1)*(fs_08/N_08) - fs_08/2;
f_16 = (0:N_16-1)*(fs_16/N_16) - fs_16/2;

figure(1);
subplot(2,1,1)
plot(f_08(1:N_08),(FT_08(1:N_08)));
grid on ;
title('Magnitude Spectrum Plot - 8kHz Sampled Speech');
xlabel('Frequency (Hz)');
ylabel('Magnitude (|H(jw)|)');

subplot(2,1,2)
plot(f_16(1:N_16),(FT_16(1:N_16)));
grid on ;
title('Magnitude Spectrum Plot - 16kHz Sampled Music');
xlabel('Frequency (Hz)');
ylabel('Magnitude (|H(jw)|)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 1 : Downsampling (M = 2)

y_08_ds = downsample(y_08,ds_factor) ; fs_08_ds = fs_08/ds_factor ;
y_16_ds = downsample(y_16,ds_factor) ; fs_16_ds = fs_16/ds_factor ;

FT_08_ds = abs(fftshift(fft(y_08_ds,length(y_08_ds)))) ;
FT_16_ds = abs(fftshift(fft(y_16_ds,length(y_16_ds)))) ;
f_08_ds = (0:length(y_08_ds)-1)*(fs_08_ds/length(y_08_ds)) - fs_08_ds/2;
f_16_ds = (0:length(y_16_ds)-1)*(fs_16_ds/length(y_16_ds)) - fs_16_ds/2;

audiowrite('DS_speech8khz.wav',y_08_ds,fs_08_ds) ;
audiowrite('DS_music16khz.wav',y_16_ds,fs_16_ds) ;

figure(2);
subplot(2,1,1) ;
plot(f_08_ds(1:length(y_08_ds)),(FT_08_ds(1:length(y_08_ds))));
grid on ;
title('Magnitude Spectrum Plot - Downsampled (M=2) Speech');
xlabel('Frequency (Hz)');
ylabel('Magnitude (|H(jw)|)');

subplot(2,1,2);
plot(f_16_ds(1:length(y_16_ds)),(FT_16_ds(1:length(y_16_ds))));
grid on ;
title('Magnitude Spectrum Plot - Downsampled (M=2) Music');
xlabel('Frequency (Hz)');
ylabel('Magnitude (|H(jw)|)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 2,3 : Equiripple Filter Design (Anti-Aliasing Filter) 
%  CAUTION : For downsampling with M = 4, change to ds_factor = 4, wp =0.22*pi, ws = 0.28*pi

y_08_lpf = filter(filter_08,y_08) ;
y_16_lpf = filter(filter_16,y_16) ;

y_08_lpf_ds = downsample(y_08_lpf,ds_factor) ; fs_08_ds = fs_08/ds_factor ;
y_16_lpf_ds = downsample(y_16_lpf,ds_factor) ; fs_16_ds = fs_16/ds_factor ;

FT_08_lpf_ds = abs(fftshift(fft(y_08_lpf_ds,length(y_08_lpf_ds)))) ;
FT_16_lpf_ds = abs(fftshift(fft(y_16_lpf_ds,length(y_16_lpf_ds)))) ;
f_08_lpf_ds = (0:length(y_08_lpf_ds)-1)*(fs_08_ds/length(y_08_lpf_ds)) - fs_08_ds/2;
f_16_lpf_ds = (0:length(y_16_lpf_ds)-1)*(fs_16_ds/length(y_16_lpf_ds)) - fs_16_ds/2;
    
audiowrite('LPF_DS_speech8khz_Q2.wav',y_08_lpf_ds,fs_08_ds) ;
audiowrite('LPF_DS_music16khz_Q2.wav',y_16_lpf_ds,fs_16_ds) ;

figure(3) ;
subplot(2,1,1) ;
plot(f_08_lpf_ds(1:length(y_08_lpf_ds)),(FT_08_lpf_ds(1:length(y_08_lpf_ds))));
grid on ;
title('Magnitude Spectrum Plot - LPF-Downsampled (M=2) Speech');
xlabel('Frequency (Hz)');
ylabel('Magnitude (|H(jw)|)');

subplot(2,1,2) ;
plot(f_16_lpf_ds(1:length(y_16_lpf_ds)),(FT_16_lpf_ds(1:length(y_16_lpf_ds))));
grid on ;
title('Magnitude Spectrum Plot - LPF-Downsampled (M=2) Music');
xlabel('Frequency (Hz)');
% ylabel('Magnitude (|H(jw)|)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 4 : Upsampling, Anti-aliasing filtering, Downsampling
%  CAUTION : For downsampling with M = 4, change to ds_factor = 4, wp =0.22*pi, ws = 0.28*pi

y_08_us = upsample(y_08,us_factor) ; fs_08_us = fs_08*us_factor ;
y_16_us = upsample(y_16,us_factor) ; fs_16_us = fs_16*us_factor ;

y_08_us_lpf = filter(filter_08,y_08_us) ;
y_16_us_lpf = filter(filter_16,y_16_us) ;

y_08_rs = downsample(y_08_us_lpf,ds_factor) ; fs_08_rs = fs_08_us/ds_factor ;
y_16_rs = downsample(y_16_us_lpf,ds_factor) ; fs_16_rs = fs_16_us/ds_factor ;

% Analysing the Upsampled Signal 

FT_08_us = abs(fftshift(fft(y_08_us,length(y_08_us)))) ;
FT_16_us = abs(fftshift(fft(y_16_us,length(y_16_us)))) ;
f_08_us = (0:length(y_08_us)-1)*(fs_08_us/length(y_08_us)) - fs_08_us/2;
f_16_us = (0:length(y_16_us)-1)*(fs_16_us/length(y_16_us)) - fs_16_us/2;

audiowrite('US_speech8khz.wav',y_08_us,fs_08_us) ;
audiowrite('US_music16khz.wav',y_16_us,fs_16_us) ;

figure(4);

subplot(2,1,1) ;
plot(f_08_us(1:length(y_08_us)),(FT_08_us(1:length(y_08_us))));
grid on ;
title('Magnitude Spectrum Plot - Upsampled (L=3) Speech');
xlabel('Frequency (Hz)');
ylabel('Magnitude (|H(jw)|)');

subplot(2,1,2) ;
plot(f_16_us(1:length(y_16_us)),(FT_16_us(1:length(y_16_us))));
grid on ;
title('Magnitude Spectrum Plot - Upsampled (L=3) Music');
xlabel('Frequency (Hz)');
ylabel('Magnitude (|H(jw)|)');

% % Analysing the Resampled Signal
 
FT_08_rs = abs(fftshift(fft(y_08_rs,length(y_08_rs)))) ;
FT_16_rs = abs(fftshift(fft(y_16_rs,length(y_16_rs)))) ;
f_08_rs = (0:length(y_08_rs)-1)*(fs_08_rs/length(y_08_rs)) - fs_08_rs/2;
f_16_rs = (0:length(y_16_rs)-1)*(fs_16_rs/length(y_16_rs)) - fs_16_rs/2;

audiowrite('RS_speech8khz.wav',y_08_rs,fs_08_rs) ;
audiowrite('RS_music16khz.wav',y_16_rs,fs_16_rs) ;

figure(5);
subplot(2,1,1) ;
plot(f_08_rs(1:length(y_08_rs)),(FT_08_rs(1:length(y_08_rs))));
grid on ;
title('Magnitude Spectrum Plot - Resampled Speech');
xlabel('Frequency (Hz)');
ylabel('Magnitude (|H(jw)|)');

subplot(2,1,2);
plot(f_16_rs(1:length(y_16_rs)),(FT_16_rs(1:length(y_16_rs))));
grid on ;
title('Magnitude Spectrum Plot - Resampled Music');
xlabel('Frequency (Hz)');
ylabel('Magnitude (|H(jw)|)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 5 : Anti-aliasing filtering -> Downsampling -> Upsampling -> Interpolation
%  CAUTION : For downsampling with M = 4, change to ds_factor = 4, wp =0.22*pi, ws = 0.28*pi

y_08_lpf = filter(filter_08,y_08) ;
y_16_lpf = filter(filter_16,y_16) ;

y_08_lpf_ds = downsample(y_08_lpf,ds_factor) ; fs_08_ds = fs_08/ds_factor ;
y_16_lpf_ds = downsample(y_16_lpf,ds_factor) ; fs_16_ds = fs_16/ds_factor ;

y_08_rs_ = upsample(y_08_lpf_ds,us_factor) ; fs_08_rs_ = fs_08_ds*us_factor ;
y_16_rs_ = upsample(y_16_lpf_ds,us_factor) ; fs_16_rs_ = fs_16_ds*us_factor ;

y_08_int = filter(interpolator,1,y_08_rs_) ;
y_16_int = filter(interpolator,1,y_16_rs_) ;

% Analysing before interpolation 

FT_08_rs_ = abs(fftshift(fft(y_08_rs_,length(y_08_rs_)))) ;
FT_16_rs_ = abs(fftshift(fft(y_16_rs_,length(y_16_rs_)))) ;
f_08_rs_ = (0:length(y_08_rs_)-1)*(fs_08_rs_/length(y_08_rs_)) - fs_08_rs_/2;
f_16_rs_ = (0:length(y_16_rs_)-1)*(fs_16_rs_/length(y_16_rs_)) - fs_16_rs_/2;

figure(6);
subplot(2,1,1) ;
plot(f_08_rs_(1:length(y_08_rs_)),(FT_08_rs_(1:length(y_08_rs_))));
grid on ;
title('Magnitude Spectrum Plot - Before Interpolation : Speech');
xlabel('Frequency (Hz)');
ylabel('Magnitude (|H(jw)|)');

subplot(2,1,2) ;
plot(f_16_rs_(1:length(y_16_rs_)),(FT_16_rs_(1:length(y_16_rs_))));
grid on ;
title('Magnitude Spectrum Plot - Before Interpolation : Music');
xlabel('Frequency (Hz)');
ylabel('Magnitude (|H(jw)|)');

% Analysing after interpolation 

audiowrite('INT_speech8khz.wav',y_08_int,fs_08_rs_) ;
audiowrite('INT_music16khz.wav',y_16_int,fs_16_rs_) ;

FT_08_int = abs(fftshift(fft(y_08_int,length(y_08_int)))) ;
FT_16_int = abs(fftshift(fft(y_16_int,length(y_16_int)))) ;
f_08_int = (0:length(y_08_int)-1)*(fs_08_rs_/length(y_08_int)) - fs_08_rs_/2;
f_16_int = (0:length(y_16_int)-1)*(fs_16_rs_/length(y_16_int)) - fs_16_rs_/2;

figure(7);
subplot(2,1,1) ;

plot(f_08_int(1:length(y_08_int)),(FT_08_int(1:length(y_08_int))));
grid on ;
title('Magnitude Spectrum Plot - After Interpolation : Speech');
xlabel('Frequency (Hz)');
ylabel('Magnitude (|H(jw)|)');

subplot(2,1,2) ;

plot(f_16_int(1:length(y_16_int)),(FT_16_int(1:length(y_16_int))));
grid on ;
title('Magnitude Spectrum Plot - After Interpolation : Music');
xlabel('Frequency (Hz)');
ylabel('Magnitude (|H(jw)|)');

function [filter]= equiripple_filter(fs,wp,ws,del,attenuation)

fp = fs*wp/(2*pi);
fst = fs*ws/(2*pi);

Hd = fdesign.lowpass('Fp,Fst,Ap,Ast',fp,fst,del,attenuation,fs);
filter = design(Hd,'equiripple');

end
