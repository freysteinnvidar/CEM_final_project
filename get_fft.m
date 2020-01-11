function [EH_fourier,s_fft] = get_fft(EH,freq)
global dt
[t_fft,s_fft] = size(EH);
Time = (0:dt:t_fft*dt - dt);
for x = 1:s_fft
    EH_fft(x,:) = fft(EH(:,x));
end
f_fft = (0:t_fft-1).*1/dt./t_fft;
[~,idx]=min(abs(f_fft-freq));
EH_fourier = EH_fft(:,idx);
