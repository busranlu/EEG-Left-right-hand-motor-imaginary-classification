%% hasta verisinin yüklenmesi
%veriyi elle ekliyoruz.
fs=200; %örnekleme frekansı

%% ön işleme
%verinin ve BBCI Toolbox'ın tanıtılması
startup_bbci_toolbox('DataDir','D:\EegMyDataDir','TmpDir','D:\TemDir');
BTB.History = 0;

% parametreler
subdir_list = {'VP001'}; %katılımcı
basename_list = {'motor_imagery1','mental_arithmetic1','motor_imagery2','mental_arithmetic1','motor_imagery3','mental_arithmetic1'}; % task type: motor imagery / recording session: 1 - 3
stimDef.eeg = {16,32; 'condition1','condition2'};

% göz gürültüsü olmadan verinin yüklenmesi
loadDir = fullfile('C:\Users\Busra\Documents\MATLAB\EEG\EEG01\subject 01\with occular artifact');
cd(loadDir);
load cnt; load mrk, load mnt; 
%sürekli eeg sinyalinin yüklenmesi


% for motor imagery: imag
cnt_temp = cnt; mrk_temp = mrk; % save data temporarily
clear cnt mrk;
[cnt.imag, mrk.imag] = proc_appendCnt({cnt_temp{1}, cnt_temp{3}, cnt_temp{5}}, {mrk_temp{1}, mrk_temp{3}, mrk_temp{5}}); % merged motor imagery cnts

% Select EEG channels only (excluding EOG channels) for classification
% clab = {'F7','FAF5','F3','AFp1','AFp2','FAF6','F4','F8','FAF1','FAF2','Cz','Pz','CFC5','CFC3','CCP5','CCP3','T7','P7','P3','PPO1','OPO1','OPO2','PPO2','P4','CFC4','CFC6','CCP4','CCP6','P8','T8','VEOG','HEOG'}
cnt.imag = proc_selectChannels(cnt.imag,'not','*EOG'); % remove EOG channels (VEOG, HEOG)
mnt.imag = mnt_setElectrodePositions(cnt.imag.clab); % update montage

% common average reference
cnt.imag = proc_commonAverageReference(cnt.imag);

% frequency band selection for common spatial pattern (CSP)
MotorChannel = {'CFC5','CFC6','CFC3','CFC4','Cz,','CCP5','CCP6','CCP3','CCP4'};
ParientalChannel = {'Pz','P3','P4','P7','P8'};
FrontalChannel = {'F7','FAF5','F3','AFp1','FAF1','AFp2','FAF2','FAF6','F4','F8'};
OccipitalChannel = {'PPO1','OPO1','OPO2','PPO2'};


%% Cheby2 bandpass filter with a passband of band_csp, with at most Rp dB of passband ripple and at least Rs dB attenuation in the stopbands that are 3 Hz wide on both sides of the passband
%0.5-50Hz 4.derece passband chby2 ile filtreleme
Wp.imag = [0.5 30]/fs/2;
Ws.imag = [0.4 35]/fs/2;
Rp.imag = 3; % in dB
Rs.imag = 30; % in dB
n=4;
[ord.imag, Ws.imag] = cheb2ord(Wp.imag, Ws.imag, Rp.imag, Rs.imag);
[filt_b.imag,filt_a.imag] = cheby2(ord.imag, Rs.imag, Ws.imag);
[filt_b.imag,filt_a.imag] = cheby2(n, Rs.imag, Ws.imag);

epo.imag = proc_filtfilt(cnt.imag, filt_b.imag, filt_a.imag);
plot(epo.imag.x)

% plot(cnt_artifact{1, 5}.x(:,1))  
%channels=cnt.imag.clab;  
% markers(:,1)=mrk.imag.event.desc;
% markers(:,2)=transpose(int64(mrk.imag.time/1000*200)); 

EEG01=epo.imag.x;
plot(EEG01(1:6000,1))

%% segmentlere ayırma
%mrks for motor imagery. 
% marker 16 = left hand imagery, 
% marker 32 = right hand imagery
markers(:,1)=mrk.imag.event.desc;
markers(:,2)=transpose(mrk.imag.time); %msec halinde örnek sayısı olmalı
for i=1:length(markers)
        markers(i,2)=int64(markers(i,2)/1000*200);
end

a=1;
for i=1:60
   if markers(i,1)==16
       markers_onalti(a,1)=markers(i,2);
       markers_onalti(a,2)=markers(i+1,2);
       a=a+1;
   end
end

b=1;
for i=1:59
   if markers(i,1)==32
       markers_otuziki(b,1)=markers(i,2);
       markers_otuziki(b,2)=markers(i+1,2);
       b=b+1;
   end
end

%sadece 1 kanal için sol el segmentlerin çıkarılması
for i=1:30
    k=1;
    for j=markers_onalti(i,1):markers_onalti(i,2)
        left_hand_imagery(i,1)=EEG01(j,27);
        k=k+1;
    end
end


%sadece 1 kanal için sağ el segmentlerin çıkarılması
for i=1:30
    l=1;
    for j=markers_otuziki(i,1):markers_otuziki(i,2)
        right_hand_imagery(i,1)=EEG01(j,27);
        l=l+1;
    end
end

%Segmentler eşit boyda olmadığı için 5490 örneklik ortak kısımlar alınacak
for i=1:30
    n=1;
    for j=1:5490
        right_hand_imagery1(i,n)=right_hand_imagery(i,j);
        left_hand_imagery1(i,n) =left_hand_imagery(i,j);
        n=n+1;
    end
end


%% feature çıkarımı
% zaman düzlemi: mean, median,variance, standard deviation, skewness, kurtosis
%ortalama
for i=1:30
    ortalama_16(i,1)=mean(left_hand_imagery1(i,:));
end

for i=1:30
    ortalama_32(i,1)=mean(right_hand_imagery1(i,:));
end

%median
for i=1:30
    medyan_16(i,1)=median(left_hand_imagery1(i,:));
end

for i=1:30
    medyan_32(i,1)=median(right_hand_imagery1(i,:));
end

%varyans
for i=1:30
    varyans_16(i,1)=var(left_hand_imagery1(i,:));
end

for i=1:30
    varyans_32(i,1)=var(right_hand_imagery1(i,:));
end

%standart sapma
for i=1:30
    std_16(i,1)=std(left_hand_imagery1(i,:));
end

for i=1:30
    std_32(i,1)=std(right_hand_imagery1(i,:));
end

%skewness skewness
for i=1:30
    ske_16(i,1)= skewness(left_hand_imagery1(i,:));
end

for i=1:30
    ske_32(i,1)= skewness(right_hand_imagery1(i,:));
end

%kurtosis 
for i=1:30
    kur_16(i,1)= kurtosis(left_hand_imagery1(i,:));
end

for i=1:30
    kur_32(i,1)= kurtosis(right_hand_imagery1(i,:));
end

%peak amplitude
for i=1:30
    peak_16(i,1)= max(left_hand_imagery1(i,:));
end

for i=1:30
    peak_32(i,1)= max(right_hand_imagery1(i,:));
end


% frekans düzlemi:
% frekans domain IRP değilse ve osilasyon varsa time domain istatiksel
% feature almıyoruz. Cünkü sinyal non gaussian
% welch ile power spectrum al sonra normal dağılım old. için istatikselleri al
% power spectrum ortalamalarından anlamlı bandı ara sonra hepsinde o bant arasındaki sinyalde featurelar cıkar

%% power spectrum
%30 segment eeg için her segmentin güç spektrumu alalım
%sol el güç spektrumu
L = 1000;        %pencere uzunluğu
noverlap = 50;   %overlap
nfft=1024;       %number of fft
fs=200;          %örnekleme frekansı
[pxx1,f1] = pwelch(left_hand_imagery1(1,:),hamming(L),noverlap,nfft,fs);
[pxx2,f2] = pwelch(left_hand_imagery1(2,:),hamming(L),noverlap,nfft,fs);
[pxx3,f3] = pwelch(left_hand_imagery1(3,:),hamming(L),noverlap,nfft,fs);
[pxx4,f4] = pwelch(left_hand_imagery1(4,:),hamming(L),noverlap,nfft,fs);
[pxx5,f5] = pwelch(left_hand_imagery1(5,:),hamming(L),noverlap,nfft,fs);
[pxx6,f6] = pwelch(left_hand_imagery1(6,:),hamming(L),noverlap,nfft,fs);

%ekrana çizdirme, sol el
figure;
subplot(3,2,1)
plot(f1(f1<10),pxx1(f1<10),'b','LineWidth',3)
xlabel('Frequency (Hz)')
ylabel('Power')
title(['LH1, Welch method, Hann. Win. Len:' num2str(L)]);
subplot(3,2,2)
plot(f2(f2<10),pxx2(f2<10),'b','LineWidth',3)
xlabel('Frequency (Hz)')
ylabel('Power')
title(['LH2, Welch method, Hann. Win. Len:' num2str(L)]);
subplot(3,2,3)
plot(f3(f3<10),pxx3(f3<10),'b','LineWidth',3)
xlabel('Frequency (Hz)')
ylabel('Power')
title(['LH3, Welch method, Hann. Win. Len:' num2str(L)]);
subplot(3,2,4)
plot(f4(f4<10),pxx4(f4<10),'b','LineWidth',3)
xlabel('Frequency (Hz)')
ylabel('Power')
title(['LH4, Welch method, Hann. Win. Len:' num2str(L)]);
subplot(3,2,5)
plot(f5(f5<10),pxx5(f5<10),'b','LineWidth',3)
xlabel('Frequency (Hz)')
ylabel('Power')
title(['LH5, Welch method, Hann. Win. Len:' num2str(L)]);
subplot(3,2,6)
plot(f6(f6<10),pxx6(f6<10),'b','LineWidth',3)
xlabel('Frequency (Hz)')
ylabel('Power')
title(['LH6, Welch method, Hann. Win. Len:' num2str(L)]);
hold on

%ekrana çizdirme, sağ el
figure;
subplot(3,2,1)
plot(f11(f11<10),pxx11(f11<10),'b','LineWidth',3)
xlabel('Frequency (Hz)')
ylabel('Power')
title(['RH1, Welch method, Hann. Win. Len:' num2str(L)]);
subplot(3,2,2)
plot(f22(f22<10),pxx22(f22<10),'b','LineWidth',3)
xlabel('Frequency (Hz)')
ylabel('Power')
title(['RH2, Welch method, Hann. Win. Len:' num2str(L)]);
subplot(3,2,3)
plot(f33(f33<10),pxx33(f33<10),'b','LineWidth',3)
xlabel('Frequency (Hz)')
ylabel('Power')
title(['RH3, Welch method, Hann. Win. Len:' num2str(L)]);
subplot(3,2,4)
plot(f44(f44<10),pxx44(f44<10),'b','LineWidth',3)
xlabel('Frequency (Hz)')
ylabel('Power')
title(['RH4, Welch method, Hann. Win. Len:' num2str(L)]);
subplot(3,2,5)
plot(f55(f55<10),pxx55(f55<10),'b','LineWidth',3)
xlabel('Frequency (Hz)')
ylabel('Power')
title(['RH5, Welch method, Hann. Win. Len:' num2str(L)]);
subplot(3,2,6)
plot(f66(f66<10),pxx66(f66<10),'b','LineWidth',3)
xlabel('Frequency (Hz)')
ylabel('Power')
title(['RH6, Welch method, Hann. Win. Len:' num2str(L)]);


%sol el, ortalama güç spektrumu
pxx_welch_sum_LH=0;
f_welch_sum_LH=0;
for i=1:30
    pxx_welch_sum_LH = pxx_welch_sum_LH + pxx1+pxx2+pxx3+pxx4+pxx5;
    f_welch_sum_LH   = f_welch_sum_LH + f1+f2+f3+f4+f5;
end
pxx_welch_mean_LH = pxx_welch_sum_LH/5;
f_welch_mean_LH   = f_welch_sum_LH/5;

plot(f_welch_mean_LH(f_welch_mean_LH<30),pxx_welch_mean_LH(f_welch_mean_LH<30),'b','LineWidth',3)
xlabel('Frequency (Hz)')
ylabel('Power')
title(['LH, Welch method mean']);


%sağ el güç spektrumu
[pxx11,f11] = pwelch(right_hand_imagery1(1,:),hamming(L),noverlap,nfft,fs);
[pxx22,f22] = pwelch(right_hand_imagery1(2,:),hamming(L),noverlap,nfft,fs);
[pxx33,f33] = pwelch(right_hand_imagery1(3,:),hamming(L),noverlap,nfft,fs);
[pxx44,f44] = pwelch(right_hand_imagery1(4,:),hamming(L),noverlap,nfft,fs);
[pxx55,f55] = pwelch(right_hand_imagery1(5,:),hamming(L),noverlap,nfft,fs);
[pxx66,f66] = pwelch(right_hand_imagery1(6,:),hamming(L),noverlap,nfft,fs);


%sağ el, ortalama güç spektrumu
for i=1:30
    pxx_welch_sum_RH = pxx_welch_sum_RH + pxxri;
    f_welch_sum_RH   = f_welch_sum_RH + fri;
end
pxx_welch_mean_RH = pxx_welch_sum_RH/30;
f_welch_mean_RH   = f_welch_sum_RH/30;

plot(f_welch_mean_RH(f_welch_mean_RH<30),pxx_welch_mean_RH(f_welch_mean_RH<30),'b','LineWidth',3)
xlabel('Frequency (Hz)')
ylabel('Power')
title(['RH, Welch method mean']);
%çıkan frekans değerleri feature kabul edilmiştir


%(frekans bantları:
% Delta δ (0.5-4),
% Theta θ (4-8), 
% Alpha α (8-14), 
% Beta β  (14-30), 
% Gamma γ (30+),
plot(f_welch_mean_RH(f_welch_mean_RH<30),pxx_welch_mean_RH(f_welch_mean_RH<30),'b','LineWidth',3)
xlim([8 14]);
xlabel('Frequency (Hz)')
ylabel('Power')
title(['EO, Welch method mean']);

 