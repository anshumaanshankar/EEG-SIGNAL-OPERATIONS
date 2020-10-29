clc
clear all
close all
f0=150; % cutoff frequency
Q=35; % Q factor 
[y, t, freqint] = plotATM('Subject00_1_edfm'); % reading the input EEG, the signal is stored in y
Fs = freqint(1); % sampling frequency that is stored in freqint is moved to Fs
figure()
subplot(3,1,1) 
plot(t,y,'b') % plotting the input signal with respect to time 
xlabel('time(s)')
ylabel('amplitude(uV)')
title('EEG pattern')

% QUANTISATION 
outmin=min(y) % calculating the minimum value of the signal 
outmax=max(y) % calculating the maximum value of the signal
num=500 % quantization levels
yquan=y/outmax(1)
delta=(outmax-outmin)/num % step size calculation
delta=delta(1)
q=delta.*[0:num-1]
q=q-((num-1)/2)*delta

subplot(3,1,2)
stem(q) % displaying the quantization levels 
xlabel('time(s)')
ylabel('Amplitude(uV)')
title('Quantisation Levels')
% loop for quantizing the wave by rounding off values to nearest quantization value
for x=1:num
yquan(find((q(x)-delta/2<=yquan)&(yquan<=q(x)+delta/2)))=q(x).*ones(1,length(find((q(x)-delta/2<=yquan)&(yquan<=q(x)+delta/2)))) 
bquan(find(yquan==q(x)))=(x-1).*ones(1,length(find(yquan==q(x))))
end
subplot(3,1,3)
plot(t,yquan, 'b') % displaying the quantized output wave
xlabel('time(s)')
ylabel('amplitude(uV)')
title('quantised waveform')
figure()
subplot(2,1,1)
plot(t,y,'r', 'linewidth',1.5) % plotting zoomed input signal
axis([0.2 0.8 -50 50]) % setting axis values to show zoomed region
title('Original EEG pattern')
xlabel('Time(s)')
ylabel('Amplitude (uV)')

subplot(2,1,2)
plot(t,yquan,'r','linewidth',1.5) % Corresponding Zoomed output to show the quantised output clearly
axis([0.2 0.8 -1 1]) % Axis setting in order to display zoomed output region
title('Quantized EEG pattern')
xlabel('Time(s)')
ylabel('Amplitude (uV)')

% DIGITISATION OF THE QUANTISED SIGNAL
n1=8;  % number of bits per sample
partition=outmin:delta:outmax;			% quantisation lines definition
codebook=outmin-(delta/2):delta:outmax+(delta/2);    % number of quantisation levels
[index,quants]=quantiz(yquan,partition,codebook);
l1=length(index); % decreasing the number of indices for one, for digitisation
convert=de2bi(index,'left-msb') 	% conversion from decimal to binary in order to perform digitisation
k=1;
for i=1:l1  % for loop showing conversion from column vector to row vector for all binary values
    for j=1:n1
        converted(k)=convert(i,j);
        j=j+1;
        k=k+1;
    end
    i=i+1;
end
figure(3);
hold on;
stairs(converted); % plots the digital version of the EEG signal
axis([0 100 -1 1]) % setting values on axes for the plot
title('DIGITAL EEG SIGNAL');
xlabel('time(s)');
ylabel('amplitude(uV)');
hold off;

% RECONSTRUCTION OF QUANTISED SIGNAL FROM THE CONVERTED SIGNAL
back=reshape(converted,n1,(length(converted)/n1));
index1=bi2de(back,'left-msb'); % conversion of binary values back into decimal
resignal=delta*index+outmin+(delta/2); % reconstructing the original signal using the quantised values

figure(4);
hold on;
plot(t,resignal); 	% plots the reconstructed signal with respect to time
axis([0.2 0.8 -1 1])	% axis setting shows that the reconstructed signal is the same as the original zoomed in output of the quantised signal
title('DEMODULATED SIGNAL');
xlabel('Time(s)');
ylabel('Amplitude (uV)');
hold off;
