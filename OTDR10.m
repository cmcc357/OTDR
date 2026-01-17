%% 空芯光纤LFM-OTDR直接检测系统仿真 (200km版本 - 优化版)
% 方案：LFM强度调制 + 直接检测 + 脉冲压缩 + 调制曲线非线性
clear; close all; clc;

fprintf('=== 200km空芯光纤LFM-OTDR直接检测仿真系统 (优化版) ===\n');

%% ====================== 系统参数设置 ======================
fprintf('设置系统参数...\n');

% 基本物理参数
c = 3e8;                    % 光速 (m/s)
n_air = 1.0;                % 空芯光纤有效折射率 (空气近似)
v = c/n_air;                % 光在光纤中速度

% LFM信号参数 (按您的要求)
PulseWidth = 20e-6;          % LFM脉冲宽度 20μs
Bandwidth = 10e6;            % LFM带宽 100MHz
% PulseWidth_TRA = PulseWidth; 
PulseWidth_TRA = 1 / Bandwidth;     % 传统脉冲宽度为LFM脉冲压缩后的宽度
f_center = 0;                % 基带信号，中心频率为0

% 系统采样参数
Fs = 20 * Bandwidth;         % 采样频率 800MHz
Ts = 1 / Fs;                   % 采样间隔 1.25ns

% 光纤参数
FiberLength = 200000;       % 光纤长度 200km
alpha_fiber = 0.2;          % 光纤衰减系数 dB/km 
alpha_linear = alpha_fiber / (10 * log10(exp(1))) / 1000; % 转换为线性系数

% 事件定义: [位置(m), 反射系数(dB), 类型(1:反射, 2:散射损耗), 宽度(m)]
Events = [...
    20000,  -60,   1,  2;      % 连接器反射
% %     45000,  -4,    2,  25;     % 轻微弯曲散射
    20500,  -65,   1,  2];     % 弱反射点
% %     60000,  -10,   2,  40];    % 损耗区域

%% ====================== 1. LFM信号生成 ======================
fprintf('1. 生成LFM探测信号...\n');

t_pulse = (0:Ts:PulseWidth-Ts)';          % LFM脉冲时间轴
t_pulse1 = (0:Ts:PulseWidth_TRA-Ts)';          % 传统脉冲时间轴
N_pulse = length(t_pulse);
mu = Bandwidth / PulseWidth;              % 调频斜率

% 生成实基带LFM信号 (强度调制信号)
LFM_Baseband = cos(pi * mu * t_pulse.^2); % 实数LFM信号
% LFM_Baseband = exp(1j * pi * mu * t_pulse.^2); % 复数LFM信号

fprintf('   LFM脉冲参数: 宽度=%.1fμs, 带宽=%.1fMHz, 时间带宽积=%.0f\n', ...
        PulseWidth*1e6, Bandwidth/1e6, Bandwidth*PulseWidth);

%% ====================== 2. MZ光强度调制仿真 (修正的非线性模型) ======================
fprintf('2. 模拟MZ光强度调制 (修正非线性模型)...\n');

% MZ调制器参数
V_pi = 4.0;                 % 半波电压 (V)
V_bias = V_pi/2;            % 偏置在正交点
Extinction_Ratio = 20;      % 消光比 20dB
ER_linear = 10^(Extinction_Ratio/10); % 线性消光比

% 驱动信号 (考虑适当的幅度)
Drive_Amplitude = 1.2;      % 驱动信号幅度 (V) - 适当减小以避免过度非线性
LFM_Drive = Drive_Amplitude * LFM_Baseband;

% 修正的MZ调制器传输函数 - 考虑调制曲线非线性
% Optical_Carrier_Power = 10e-3;            % 激光器输出光功率 10mW

% MZ调制器理想传输函数: P_out/P_in = 0.5 * [1 + cos(π*(V_bias + V_drive)/V_pi)]
% 添加非线性失真：实际调制器响应与理想余弦函数的偏差

Modulator_Input = V_bias + LFM_Drive;

% 理想余弦响应
ideal_response = 0.5 * (1 + cos(pi * Modulator_Input / V_pi));

% 添加调制曲线非线性失真 (模拟实际MZ调制器的非线性)
% 使用多项式拟合实际调制器响应与理想响应的偏差
nonlinear_coeff = [0.02, -0.01, 0.005]; % 非线性系数
nonlinear_distortion = nonlinear_coeff(1) * Modulator_Input.^2 + ...
                      nonlinear_coeff(2) * Modulator_Input.^3 + ...
                      nonlinear_coeff(3) * Modulator_Input.^4;

% 实际调制器响应 (理想响应 + 非线性失真)
actual_response = ideal_response + nonlinear_distortion;

% 频域相关实现匹配滤波
N_corr = length(Modulator_Input) + length(Modulator_Input) - 1;
S_rx = fft(LFM_Drive, N_corr);
S_ref = fft(LFM_Drive, N_corr);
S_corr = S_rx .* conj(S_ref);
Modulator_Input_compress = ifft(S_corr);

% % 时域匹配
% matched_filter = conj(fliplr(actual_response));
% Modulator_Input_compress = conv(actual_response, matched_filter);
Modulator_Input_compress = Modulator_Input_compress(1:length(Modulator_Input));
t_receive1 = (0:Ts:(N_corr-1)*Ts);
figure;
subplot(2,1,1);
plot(t_pulse*1e6, LFM_Baseband, 'b-', 'LineWidth', 1.5);
xlabel('时间 (μs)'); ylabel('幅度 (V)');
title('LFM信号 (时域)');
grid on;
subplot(2,1,2);
plot(t_pulse*1e6, abs(Modulator_Input_compress), 'g-','LineWidth', 1.5);
title('LFM自相关信号(时域)')
xlabel('时间 (us)'); ylabel('幅度(V');
grid on;

% 频域相关实现匹配滤波
N_corr = length(Modulator_Input) + length(Modulator_Input) - 1;
actual_response_test = actual_response;
actual_response_test = actual_response_test - mean(actual_response_test);
% N_corr_1 = 2^nextpow2(N_corr);
% S_rx = fft(actual_response_test, N_corr_1);
% S_ref = fft(actual_response_test, N_corr_1);
% S_corr = S_rx .* conj(S_ref);
% Modulator_Input_compress = ifft(S_corr);

% % 时域匹配
matched_filter = conj(fliplr(actual_response_test));
matched_filter = [];
for i = 1:length(actual_response_test)
matched_filter(i) = actual_response_test(end - i+1);
end
Modulator_Input_compress = conv(actual_response_test, matched_filter);
% % Modulator_Input_compress = Modulator_Input_compress(1:length(Modulator_Input));
t_receive2 = (0:Ts:(N_corr-1)*Ts);
% t_receive3 = (0:Ts:(N_corr_1-1)*Ts);
figure;
subplot(2,1,1);
plot(t_pulse*1e6, actual_response_test, 'b-', 'LineWidth', 1.5);
xlabel('时间 (μs)'); ylabel('幅度 (V)');
title('调制光功率信号 (时域)');
grid on;
subplot(2,1,2);
plot(t_receive2*1e6, abs(Modulator_Input_compress), 'g-','LineWidth', 1.5);
title('调制光功率自相关信号(时域)')
xlabel('时间 (us)'); ylabel('幅度(V');
grid on;

% actual_response = ideal_response;
% 应用消光比和功率归一化
% min_power = 1 / ER_linear;
% actual_response = max(actual_response, min_power);
% figure;
% plot(t_pulse,actual_response);
% p_LFM = sum(actual_response.^2 / PulseWidth / Fs);
p_LFM_peak = max(abs(actual_response));
% fprintf('LFM光信号功率：%.8f\n',p_LFM);
% fprintf('LFM光信号峰值功率：%.8f\n',p_LFM_peak);
p_LFM_peak1 = 100*10^-3; % LFM光信号峰值功率为100mW
Optical_Carrier_Power = p_LFM_peak1 / p_LFM_peak;
LFM_Optical = Optical_Carrier_Power * actual_response;
p_LFM_Optical_peak = (max(abs(LFM_Optical)));
fprintf('LFM光信号峰值功率：%.8f\n',p_LFM_Optical_peak);
p_lfm = mean(LFM_Optical);
fprintf('LFM光信号功率：%.8f\n',p_lfm);
% figure;
% plot(t_pulse*1e6, LFM_Optical, 'r-', 'LineWidth', 1.5);
% xlabel('时间 (μs)'); ylabel('调制光功率 (W)');
% title('LFM调制光功率信号 (时域)');
% grid on;
% 计算调制深度和线性度
modulation_depth = (max(LFM_Optical) - min(LFM_Optical)) / Optical_Carrier_Power;
fprintf('   调制深度: %.1f%%, 驱动幅度: %.2fV\n', modulation_depth*100, Drive_Amplitude);

%% ====================== 3. 光纤信道建模 ======================
fprintf('3. 构建200km空芯光纤信道模型...\n');

% 计算接收时间窗口 (200km往返时间约1.33ms)
max_delay = 2 * FiberLength / v;          % 最大往返时延
t_receive = (0:Ts:max_delay)';            % 接收时间轴
N_receive = length(t_receive);
t_receive = (0:Ts:(N_receive-1)*Ts)';
distance_axis = t_receive * v / 2;        % 距离轴 (米)

% 初始化光纤冲激响应 (功率域)
h_fiber = zeros(N_receive, 1);

% 构建空芯光纤瑞利散射背景 (比实芯光纤弱约20dB)
Rayleigh_Level_HCF = -100;                % 200km空芯光纤瑞利散射水平 (dB)
% rayleigh_background = (10^(Rayleigh_Level_HCF/20)) * ...
%                      exp(-alpha_linear * distance_axis) .* ...
%                      (0.3 + 0.7*randn(N_receive, 1));
rayleigh_background = (10^(Rayleigh_Level_HCF/10)) * ...
                     exp(-2*alpha_linear * distance_axis) * (v * Ts / 2);           
% 模拟瑞利背向散射信号的随机波动，代表了由于光纤的不均匀性和环境变化引起的信号起伏，0.3和0.7是工程经验值，可根据不同光纤微调
% rayleigh_background = (10^(Rayleigh_Level_HCF/20)) * ...
%                      exp(-alpha_linear * distance_axis) .* ...
%                      raylrnd(0.25, N_receive, 1);
                 
h_fiber = h_fiber + rayleigh_background;

% 添加离散事件
time_idex_all = [];
for i = 1:size(Events, 1)
    pos = Events(i, 1);
    loss_dB = Events(i, 2);
    event_type = Events(i, 3);
    event_width = Events(i, 4);
    
    % 计算到事件位置的衰减
    attenuation = exp(-alpha_linear * 2 * pos);
    coeff_linear = 10^(loss_dB/10) * attenuation;
    
    % 找到事件对应的时间索引
    time_idx = round(pos / (v * Ts / 2));
    if time_idx > N_receive
        continue;
    end
    
    switch event_type
        case 1  % 反射事件 (狄拉克脉冲)
            h_fiber(time_idx) = h_fiber(time_idx) + coeff_linear;
            
        case 2  % 散射损耗事件 (高斯形状)
            width_samples = round(event_width / (v * Ts / 2));
            idx_start = max(1, time_idx - 2*width_samples);
            idx_end = min(N_receive, time_idx + 2*width_samples);
            x = (idx_start:idx_end)' - time_idx;
            gaussian_window = exp(-(x.^2) / (2*(width_samples/2)^2));
            h_fiber(idx_start:idx_end) = h_fiber(idx_start:idx_end) + ...
                                       coeff_linear * gaussian_window;
    end
    time_idex_all = [time_idex_all, time_idx];
end

%% ====================== 4. 信号传输与APD直接检测 ======================
fprintf('4. 模拟信号传输与APD直接检测...\n');

% 信号传输 (卷积过程)
% Tx_Optical_Power = LFM_Baseband;
Tx_Optical_Power = LFM_Optical;           % 发射光功率
Rx_Optical_Signal = conv(Tx_Optical_Power, h_fiber);
Rx_Optical_Signal = Rx_Optical_Signal(1:N_receive); % 裁剪
figure;
plot(t_receive*1e6, Rx_Optical_Signal, 'b-', 'LineWidth', 1.5);
xlabel('时间 (μs)'); ylabel('功率(W)');
title('光纤输出信号 (时域)');
grid on;
% Rx_Optical_Signal = Rx_Optical_Signal .* exp(-alpha_linear * distance_axis);

% APD参数 (实际参数)
APD_Responsivity = 0.9;                   % APD响应度 A/W (1550nm)
APD_Gain = 26;                            % APD增益
APD_Dark_Current = 5e-9;                  % 暗电流 5nA
APD_Excess_Noise_Factor = 13;             % 过剩噪声因子
APD_Bandwidth = Bandwidth;                    % APD带宽B
APD_Bandwidth_TRA = Bandwidth;                 % 传统APD带宽B
% 光电转换
Detected_Current = APD_Responsivity * APD_Gain * Rx_Optical_Signal;
fprintf('光电转换后的信号功率=%.2e W\n',(Detected_Current(1)));

% APD噪声计算 (实际噪声模型)
% 1. 散粒噪声 (信号相关)
Shot_Noise_Power = 2 * 1.6e-19 * (Detected_Current + APD_Dark_Current) * ...
                   APD_Gain^2 * APD_Excess_Noise_Factor * APD_Bandwidth;

% 2. 热噪声 (与信号无关)
Temperature = 300;                        % 温度 300K
Load_Resistance = 50;                     % 负载电阻 50Ω
Thermal_Noise_Power = 4 * 1.38e-23 * Temperature * APD_Bandwidth / Load_Resistance;

% 总噪声标准差
% Total_Noise_Power = 2e-12;
Total_Noise_Power = Shot_Noise_Power + Thermal_Noise_Power;
Noise_Std = sqrt(Total_Noise_Power);
fprintf('平均前总噪声=%.2e W\n', mean(Total_Noise_Power));

% 平均降噪
N = 10^5; % 平均次数
Noise_Std = Noise_Std /sqrt(N);
Total_Noise_Power = Total_Noise_Power / N;

% 加入实际APD噪声
Rx_Electrical_Noisy = Detected_Current + Noise_Std .* randn(size(Detected_Current));
fprintf('加噪声后的信号功率=%.2e W\n',(Rx_Electrical_Noisy(1)));

fprintf('   APD噪声统计: 散粒噪声=%.2e A, 热噪声=%.2e A, 平均后总噪声=%.2e W\n', ...
        sqrt(mean(Shot_Noise_Power)), sqrt(mean(Thermal_Noise_Power)), mean(Total_Noise_Power));

%% ====================== 5. 脉冲压缩处理 - 两种方法比较 ======================
fprintf('5. 进行脉冲压缩处理 - 两种方法比较...\n');

% 5.1 匹配滤波方法
fprintf('  a) 匹配滤波方法...\n');
N_corr = length(Rx_Electrical_Noisy) + length(LFM_Baseband) + - 1;
N_corr1 = 2^nextpow2(N_corr);
% 去直流分量
actual_response_r = actual_response;
Rx_Electrical_Noisy1 = Rx_Electrical_Noisy;
actual_response_r = actual_response_r - mean(actual_response_r);
Rx_Electrical_Noisy = Rx_Electrical_Noisy - mean(Rx_Electrical_Noisy);
% 频域相关实现匹配滤波
S_rx = fft(Rx_Electrical_Noisy, N_corr1);
LFM_Baseband = actual_response_r;
S_ref = fft(LFM_Baseband, N_corr1);
S_corr = S_rx .* conj(S_ref);
Compressed_MF = ifft(S_corr);
% 时域匹配滤波
% matched_filter_lfm = [];
% for i =1:length(actual_response_r)
%     matched_filter_lfm(i) = actual_response_r(end - i +1);
% end
% Compressed_MF = conv(matched_filter_lfm , Rx_Electrical_Noisy);
Compressed_MF = Compressed_MF(1:length(Rx_Electrical_Noisy));
fprintf('压缩后的信号功率=%.2e W\n',(Compressed_MF(1)));
% Compressed_MF = Rx_Electrical_Noisy1;

% % 5.2 分数阶傅里叶变换方法 (FrFT)
% fprintf('  b) 分数阶傅里叶变换方法...\n');
% FrFT_available = false;
% 
% % 检查FrFT工具箱是否可用，如果不可用则使用替代实现
% if exist('frft', 'file') == 2
%     FrFT_available = true;
%     % 使用FrFT工具箱
%     p_optimal = 2 * asec(sqrt(4 + (Bandwidth/PulseWidth)^2)) / pi; % FrFT的最优分阶数
%     Compressed_FrFT = frft(Rx_Electrical_Noisy, p_optimal);
% else
%     fprintf('    警告: 未找到FrFT工具箱，使用改进的STFT方法替代\n');
%     % 使用改进的短时傅里叶变换作为替代
%     Compressed_FrFT = improved_stft_compression(Rx_Electrical_Noisy, LFM_Baseband, Fs, mu);
% end

% 包络检测和对数转换
Compressed_MF_dB = 10*log10(abs(Compressed_MF));
% Compressed_FrFT_dB = 10*log10(abs(Compressed_FrFT));

%% ====================== 6. 传统脉冲OTDR对比 ======================
fprintf('6. 生成传统脉冲OTDR对比...\n');

% 生成传统矩形脉冲 (相同能量)
RectPulse_Energy = sum(LFM_Baseband.^2);  % LFM脉冲能量
RectPulse_Amplitude = sqrt(RectPulse_Energy / PulseWidth_TRA / Fs);
RectPulse = RectPulse_Amplitude * ones(size(t_pulse1));
% 传统脉冲光调制（MZ）
% 添加非线性失真：实际调制器响应与理想余弦函数的偏差
RectPulse_Drive = Drive_Amplitude * RectPulse;
Modulator_Input1 = V_bias + RectPulse_Drive;
% figure;
% plot(t_pulse1*1e6, RectPulse_Drive, 'b-', 'LineWidth', 1.5);
% xlabel('时间 (μs)'); ylabel('驱动电压 (V)');
% title('传统驱动信号 (时域)');
% grid on;
% 理想余弦响应
ideal_response1 = 0.5 * (1 + cos(pi * Modulator_Input1 / V_pi));
% 添加调制曲线非线性失真 (模拟实际MZ调制器的非线性)
% 使用多项式拟合实际调制器响应与理想响应的偏差
nonlinear_coeff = [0.02, -0.01, 0.005]; % 非线性系数
nonlinear_distortion1 = nonlinear_coeff(1) * Modulator_Input1.^2 + ...
                      nonlinear_coeff(2) * Modulator_Input1.^3 + ...
                      nonlinear_coeff(3) * Modulator_Input1.^4;

% 实际调制器响应 (理想响应 + 非线性失真)
actual_response1 = ideal_response1 + nonlinear_distortion1;
% actual_response1 = ideal_response1;
% 应用消光比和功率归一化
% min_power = 1 / ER_linear;
% actual_response1 = max(actual_response1, min_power);
% figure;
% plot(t_pulse,actual_response1);
p_TRA_peak = (max(abs(actual_response1)));
Optical_Carrier_Power1 = p_LFM_peak1 / p_TRA_peak;
Trad_Optical = Optical_Carrier_Power1 * actual_response1;
p_TRA_peak1 = (max(abs(Trad_Optical)));
p_TRA1 = mean(Trad_Optical);
fprintf('传统光信号峰值功率：%.8f\n',p_TRA_peak1);
fprintf('传统光信号功率：%.8f\n',p_TRA1);
% figure;
% plot(t_pulse1*1e6, Trad_Optical, 'r-', 'LineWidth', 1.5);
% xlabel('时间 (μs)'); ylabel('光调制功率 (V)');
% title('传统光调制功率信号 (时域)');
% grid on;
% 传统OTDR信号处理
Traditional_Trace = conv(Trad_Optical, h_fiber);
Traditional_Trace = Traditional_Trace(1:N_receive);
Traditional_Current = APD_Responsivity * APD_Gain * Traditional_Trace;
% APD噪声计算 (实际噪声模型)
% 1. 散粒噪声 (信号相关)
Shot_Noise_Power_tra = 2 * 1.6e-19 * (Traditional_Current + APD_Dark_Current) * ...
                   APD_Gain^2 * APD_Excess_Noise_Factor * APD_Bandwidth_TRA;

% 2. 热噪声 (与信号无关)
Temperature = 300;                        % 温度 300K
Load_Resistance = 50;                     % 负载电阻 50Ω
Thermal_Noise_Power_tra = 4 * 1.38e-23 * Temperature * APD_Bandwidth_TRA / Load_Resistance;

% 总噪声标准差
% Total_Noise_Power = 2e-12;
Total_Noise_Power_tra = Shot_Noise_Power_tra + Thermal_Noise_Power_tra;
Noise_Std_tra = sqrt(Total_Noise_Power_tra);
fprintf('传统平均前总噪声=%.2e W\n', mean(Total_Noise_Power_tra));

% 平均降噪
N = 10^5; % 平均次数
Noise_Std_tra = Noise_Std_tra /sqrt(N);
Total_Noise_Power_tra = Total_Noise_Power_tra / N;
fprintf('   传统APD噪声统计: 散粒噪声=%.2e A, 热噪声=%.2e A, 传统平均后总噪声=%.2e W\n', ...
        sqrt(mean(Shot_Noise_Power_tra)), sqrt(mean(Thermal_Noise_Power_tra)), mean(Total_Noise_Power_tra));
Traditional_Noisy = Traditional_Current + Noise_Std_tra .* randn(size(Traditional_Current));
% Traditional_Noisy = Traditional_Current + Noise_Std .* randn(size(Traditional_Current));

N_corr_tra = length(Traditional_Noisy) + length(RectPulse) - 1;
% 去直流
% Traditional_Noisy = Traditional_Noisy - mean(Traditional_Noisy);
% actual_response1 = actual_response1 - mean(actual_response1);
% 频域相关实现匹配滤波
S_rx_tra = fft(Traditional_Noisy, N_corr_tra);
RectPulse = actual_response1;
S_ref_tra = fft(RectPulse, N_corr_tra);
S_corr_tra = S_rx_tra .* conj(S_ref_tra);
Compressed_tra = ifft(S_corr_tra);
% 能量归一化
% ref_energy = sum(abs(LFM_Baseband).^2);
% Compressed_MF = Compressed_MF / ref_energy;
Compressed_tra = Compressed_tra(1:length(Traditional_Noisy));
% Compressed_tra = Traditional_Noisy;
Traditional_dB = 10*log10(abs(Compressed_tra));

%% ====================== 7. 性能分析 ======================
fprintf('7. 性能分析...\n');

% 计算动态范围
a = round(15000 / (v / 2 * Ts));
d = round(3000 / (v / 2 * Ts));
signal_region_MF = Compressed_MF_dB(1:a);
signal_MF = max(signal_region_MF);
b = round(0.2 * FiberLength / (v / 2 * Ts));
% noise_region_MF = Compressed_MF_dB(end-b:end-d);
noise_region_MF = Compressed_MF_dB(end-b:end);
noise_floor_MF = mean(noise_region_MF);
dynamic_range_MF =signal_MF - noise_floor_MF;

% signal_region_FrFT = Compressed_FrFT_dB(1:a);
% signal_FrFT = max(signal_region_FrFT);
% noise_region_FrFT = Compressed_FrFT_dB(end-b:end);
% noise_floor_FrFT = mean(noise_region_FrFT);
% dynamic_range_FrFT = signal_FrFT - noise_floor_FrFT;

signal_region_TRA = Traditional_dB(1:a);
signal_TRA = max(signal_region_TRA);
noise_region_Trad = Traditional_dB(end-b:end);
noise_floor_Trad = mean(noise_region_Trad);
dynamic_range_Trad =signal_TRA - noise_floor_Trad;


theoretical_resolution = c / (2 * n_air * Bandwidth); % 理论分辨率

% 计算LFM（压缩前）实际分辨率 (基于反射事件的3dB宽度)
pulse_compressed1 = abs(Rx_Electrical_Noisy1(1: length(Rx_Electrical_Noisy1)));
pulse_compressed1 = pulse_compressed1 / max(pulse_compressed1);
half_power1 = max(pulse_compressed1)/sqrt(2);
above_half1 = pulse_compressed1 > half_power1;
num_band1 = sum(above_half1);
fprintf('    3dBLFM脉冲(未压缩)宽度采样个数: %.1f \n', num_band1);
pulse_width_s1 = num_band1 * Ts;
achieved_resolution_MF1 = pulse_width_s1 * c / 2;

% 计算LFM（压缩后）实际分辨率 (基于反射事件的3dB宽度)
pulse_compressed = abs(Compressed_MF(1: length(Compressed_MF)));
pulse_compressed = pulse_compressed / max(pulse_compressed);
half_power = max(pulse_compressed)/sqrt(2);
above_half = pulse_compressed > half_power;
num_band = sum(above_half);
fprintf('    3dBLFM脉冲(压缩后)宽度采样个数: %.1f \n', num_band);
pulse_width_s = num_band * Ts;
achieved_resolution_MF = pulse_width_s * c / 2;

% 计算传统脉冲实际分辨率 (基于反射事件的3dB宽度)
pulse_compressed_tra = abs(Compressed_tra(1:length(Compressed_tra)));
pulse_compressed_tra = pulse_compressed_tra / max(pulse_compressed_tra);
half_power_tra = max(pulse_compressed_tra)/sqrt(2);
above_half_tra = pulse_compressed_tra > half_power_tra;
num_band_tra = sum(above_half_tra);
fprintf('    3dB传统脉冲宽度采样个数: %.1f \n', num_band_tra);
pulse_width = num_band_tra * Ts;
achieved_resolution_tra = pulse_width * c / 2;

% pulse_compressed_FrFT = abs(Compressed_FrFT(1:min(2000, length(Compressed_FrFT))));
% pulse_compressed_FrFT = pulse_compressed_FrFT / max(pulse_compressed_FrFT);
% above_half_FrFT = pulse_compressed_FrFT > (max(pulse_compressed_FrFT))/sqrt(2);
% pulse_width_ns_FrFT = sum(above_half_FrFT) * Ts * 1e9;
% achieved_resolution_FrFT = pulse_width_ns_FrFT * 1e-9 * c / 2;


fprintf('  性能指标:\n');
fprintf('    理论分辨率: %.1f m\n', theoretical_resolution);
fprintf('    MF方法动态范围: %.1f dB, 实测分辨率: %.1f m\n', dynamic_range_MF, achieved_resolution_MF);
% if FrFT_available
%     fprintf('    FrFT方法动态范围: %.1f dB, 实测分辨率: %.1f m\n', dynamic_range_FrFT, achieved_resolution_FrFT);
% else
%     fprintf('    STFT方法动态范围: %.1f dB, 实测分辨率: %.1f m\n', dynamic_range_FrFT, achieved_resolution_FrFT);
% end
fprintf('    传统OTDR动态范围: %.1f dB, 实测分辨率: %.1f m\n', dynamic_range_Trad, achieved_resolution_tra);

%% ====================== 8. 结果可视化 ======================
fprintf('8. 生成结果图表...\n');

% 图1: MZ调制器非线性特性分析
figure('Position', [100, 100, 1400, 900], 'Name', 'MZ调制器非线性特性分析');

subplot(2,3,1);
plot(t_pulse*1e6, LFM_Drive, 'b-', 'LineWidth', 1.5);
xlabel('时间 (μs)'); ylabel('驱动电压 (V)');
title('LFM驱动信号 (时域)');
grid on;

subplot(2,3,2);
% 显示MZ调制器非线性特性 - 修正版本
V_test = linspace(-V_pi, V_pi, 1000);
ideal_response_test = 0.5 * (1 + cos(pi * (V_bias + V_test) / V_pi));
nonlinear_distortion_test = nonlinear_coeff(1) * (V_bias + V_test).^2 + ...
                           nonlinear_coeff(2) * (V_bias + V_test).^3 + ...
                           nonlinear_coeff(3) * (V_bias + V_test).^4;
actual_response_test = ideal_response_test + nonlinear_distortion_test;
% actual_response_test = max(actual_response_test, 1/ER_linear);

plot(V_test, ideal_response_test, 'b-', 'LineWidth', 2);
hold on;
plot(V_test, actual_response_test, 'r-', 'LineWidth', 2);
plot(LFM_Drive, actual_response, 'g-', 'LineWidth', 2);

% 获取LFM_Drive的范围
LFM_min = min(LFM_Drive);
LFM_max = max(LFM_Drive);

% 获取y轴范围用于定位文本和虚线
ylimits = ylim;
y_text = ylimits(1) + 0.05 * (ylimits(2) - ylimits(1));  % 文本位置在y轴底部5%处
y_line = ylimits(2);  % 虚线延伸到整个y轴范围

% 添加垂直虚线标记LFM_Drive范围
plot([LFM_min, LFM_min], [ylimits(1), y_line], 'k--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
plot([LFM_max, LFM_max], [ylimits(1), y_line], 'k--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);

% % 在横坐标轴上添加marker标记端点
% plot(LFM_min, ylimits(1), 'kv', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
% plot(LFM_max, ylimits(1), 'kv', 'MarkerSize', 8, 'MarkerFaceColor', 'red');

% 更新横坐标刻度以显示端点值
current_xticks = xticks;
% 确保包含LFM_Drive的端点值
new_xticks = unique([current_xticks, LFM_min, LFM_max]);
xticks(new_xticks);

xlabel('驱动电压偏移 (V)'); ylabel('归一化光功率');
title('MZ调制器传输特性');
legend('理想响应', '实际响应(含非线性)', '工作点', 'Location', 'best'); 
grid on;

subplot(2,3,3);
plot(t_pulse*1e6, LFM_Optical, 'r-', 'LineWidth', 1.5);
xlabel('时间 (μs)'); ylabel('调制光功率 (W)');
title('LFM调制光功率信号 (含非线性)');
grid on;

subplot(2,3,4);
[Pxx, f] = pwelch(LFM_Drive, [], [], [], Fs);
plot(f/1e6, 10*log10(Pxx/max(Pxx)), 'r-', 'LineWidth', 1.5);
xlabel('频率 (MHz)'); ylabel('归一化功率谱密度 (dB)');
title('LFM驱动信号频谱');
xlim([-2*Bandwidth/10^6, 2*Bandwidth/10^6]); grid on;

subplot(2,3,5);
plot(t_pulse1*1e6, RectPulse_Drive, 'b-', 'LineWidth', 1.5);
xlabel('时间 (μs)'); ylabel('驱动电压 (V)');
title('传统驱动信号 (时域)');
grid on;

subplot(2,3,6);
plot(t_pulse1*1e6, Trad_Optical, 'r-', 'LineWidth', 1.5);
xlabel('时间 (μs)'); ylabel('光调制功率 (V)');
title('传统光调制功率信号 (含非线性)');
grid on;

% 脉冲压缩结果图
figure('Position', [100, 100, 1400, 900], 'Name', 'LFM脉冲压缩结果');
subplot(2,2,1);
plot(t_receive*10^6, 10*log10(abs(Rx_Electrical_Noisy1)), 'b-','LineWidth', 1.5);
title('LFM未压缩信号')
xlabel('时间 (us)'); ylabel('幅度(A)');
grid on;
for i = 1: length(time_idex_all)
n = time_idex_all(i); 
% 2. 提取该点的 X 和 Y 坐标
% 注意：必须保持和 plot 函数中完全一致的变换
x_val = t_receive(n) * 10^6;      % X轴坐标
y_val = 10*log10(abs(Rx_Electrical_Noisy1(n)));    % Y轴坐标

% 3. 锁定当前图层，防止覆盖
hold on; 

% 4. (可选) 在该点画一个红色的圆圈标记，使其更明显
plot(x_val, y_val, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);

% 5. 添加箭头和数值文本
% str 是显示的字符串：\leftarrow 生成左箭头，num2str 将数值转为字符串
str = [' \leftarrow ' num2str(y_val)]; 
text(x_val, y_val, str, 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'red');

% 6. 解除锁定
hold off;
end
subplot(2,2,2);
[Pxx, f] = pwelch(Rx_Electrical_Noisy1, [], [], [], Fs);
plot(f/1e6, 10*log10(Pxx/max(Pxx)), 'r-', 'LineWidth', 1.5);
xlabel('频率 (MHz)'); ylabel('归一化功率谱密度 (dB)');
title('LFM未压缩信号频谱');
grid on;
subplot(2,2,3);
plot(t_receive*10^6, 10*log10(abs(Compressed_MF)), 'g-','LineWidth', 1.5);
title('LFM压缩信号')
xlabel('时间 (us)'); ylabel('幅度(A)');
grid on;
% 1. 指定你要标注的第 n 个点
for i = 1: length(time_idex_all)
n = time_idex_all(i); 
% 2. 提取该点的 X 和 Y 坐标
% 注意：必须保持和 plot 函数中完全一致的变换
x_val = t_receive(n) * 10^6;      % X轴坐标
y_val = 10*log10(abs(Compressed_MF(n)));    % Y轴坐标

% 3. 锁定当前图层，防止覆盖
hold on; 

% 4. (可选) 在该点画一个红色的圆圈标记，使其更明显
plot(x_val, y_val, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);

% 5. 添加箭头和数值文本
% str 是显示的字符串：\leftarrow 生成左箭头，num2str 将数值转为字符串
str = [' \leftarrow ' num2str(y_val)]; 
text(x_val, y_val, str, 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'red');

% 6. 解除锁定
hold off;
end
subplot(2,2,4);
[Pxx, f] = pwelch(Compressed_MF, [], [], [], Fs);
plot(f/1e6, 10*log10(Pxx/max(Pxx)), 'r-', 'LineWidth', 1.5);
xlabel('频率 (MHz)'); ylabel('归一化功率谱密度 (dB)');
title('LFM压缩信号频谱');
grid on;

figure('Position', [100, 100, 1400, 900], 'Name', '传统脉冲压缩结果');
subplot(2,2,1);
plot(t_receive*10^6, 10*log10(abs(Traditional_Noisy)), 'b-','LineWidth', 1.5);
title('传统未压缩信号')
xlabel('时间 (us)'); ylabel('幅度(A)');
grid on;
for i = 1: length(time_idex_all)
n = time_idex_all(i); 
% 2. 提取该点的 X 和 Y 坐标
% 注意：必须保持和 plot 函数中完全一致的变换
x_val = t_receive(n) * 10^6;      % X轴坐标
y_val = 10*log10(abs(Traditional_Noisy(n)));    % Y轴坐标

% 3. 锁定当前图层，防止覆盖
hold on; 

% 4. (可选) 在该点画一个红色的圆圈标记，使其更明显
plot(x_val, y_val, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);

% 5. 添加箭头和数值文本
% str 是显示的字符串：\leftarrow 生成左箭头，num2str 将数值转为字符串
str = [' \leftarrow ' num2str(y_val)]; 
text(x_val, y_val, str, 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'red');

% 6. 解除锁定
hold off;
end
subplot(2,2,2);
[Pxx, f] = pwelch(Traditional_Noisy, [], [], [], Fs);
plot(f/1e6, 10*log10(Pxx/max(Pxx)), 'r-', 'LineWidth', 1.5);
xlabel('频率 (MHz)'); ylabel('归一化功率谱密度 (dB)');
title('传统脉冲未压缩信号频谱');
grid on;
subplot(2,2,3);
plot(t_receive*10^6, 10*log10(abs(Compressed_tra)), 'g-','LineWidth', 1.5);
title('传统压缩信号')
xlabel('时间 (us)'); ylabel('幅度(A)');
grid on;
for i = 1: length(time_idex_all)
n = time_idex_all(i); 
% 2. 提取该点的 X 和 Y 坐标
% 注意：必须保持和 plot 函数中完全一致的变换
x_val = t_receive(n) * 10^6;      % X轴坐标
y_val = 10*log10(abs(Compressed_tra(n)));    % Y轴坐标

% 3. 锁定当前图层，防止覆盖
hold on; 

% 4. (可选) 在该点画一个红色的圆圈标记，使其更明显
plot(x_val, y_val, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);

% 5. 添加箭头和数值文本
% str 是显示的字符串：\leftarrow 生成左箭头，num2str 将数值转为字符串
str = [' \leftarrow ' num2str(y_val)]; 
text(x_val, y_val, str, 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'red');

% 6. 解除锁定
hold off;
end
subplot(2,2,4);
[Pxx, f] = pwelch(Compressed_tra, [], [], [], Fs);
plot(f/1e6, 10*log10(Pxx/max(Pxx)), 'r-', 'LineWidth', 1.5);
xlabel('频率 (MHz)'); ylabel('归一化功率谱密度 (dB)');
title('传统脉冲压缩信号频谱');
grid on;

% 图2: 脉冲压缩方法比较
% figure('Position', [150, 150, 1400, 1000], 'Name', '脉冲压缩方法比较分析');

% subplot(3,2,1);
% zoom_start = 0; zoom_end = 50;
% zoom_idx = distance_axis/1000 >= zoom_start & distance_axis/1000 <= zoom_end;
% plot(distance_axis(zoom_idx)/1000, Compressed_MF_dB(zoom_idx), 'b-', 'LineWidth', 1.5);
% hold on;
% plot(distance_axis(zoom_idx)/1000, Compressed_FrFT_dB(zoom_idx), 'r-', 'LineWidth', 1.5);
% xlabel('距离 (km)'); ylabel('相对强度 (dB)');
% title('脉冲压缩结果对比 (0-50km)');
% if FrFT_available
%     legend('匹配滤波', 'FrFT', 'Location', 'southwest');
% else
%     legend('匹配滤波', 'STFT替代', 'Location', 'southwest');
% end
% xlim([zoom_start, zoom_end]); ylim([-100, 5]); grid on;

% subplot(3,2,2);
% zoom_start = 45; zoom_end = 55; % 50km事件附近
% zoom_idx = distance_axis/1000 >= zoom_start & distance_axis/1000 <= zoom_end;
% plot(distance_axis(zoom_idx)/1000, Compressed_MF_dB(zoom_idx), 'b-', 'LineWidth', 2);
% hold on;
% plot(distance_axis(zoom_idx)/1000, Compressed_FrFT_dB(zoom_idx), 'r-', 'LineWidth', 1.5);
% xlabel('距离 (km)'); ylabel('相对强度 (dB)');
% title('50km事件细节');
% if FrFT_available
%     legend('匹配滤波', 'FrFT', 'Location', 'southwest');
% else
%     legend('匹配滤波', 'STFT替代', 'Location', 'southwest');
% end
% xlim([zoom_start, zoom_end]); ylim([-100, 5]); grid on;

% subplot(3,2,3);
% zoom_start = 95; zoom_end = 105; % 100km事件附近
% zoom_idx = distance_axis/1000 >= zoom_start & distance_axis/1000 <= zoom_end;
% plot(distance_axis(zoom_idx)/1000, Compressed_MF_dB(zoom_idx), 'b-', 'LineWidth', 2);
% hold on;
% plot(distance_axis(zoom_idx)/1000, Compressed_FrFT_dB(zoom_idx), 'r-', 'LineWidth', 1.5);
% xlabel('距离 (km)'); ylabel('相对强度 (dB)');
% title('100km事件细节');
% if FrFT_available
%     legend('匹配滤波', 'FrFT', 'Location', 'southwest');
% else
%     legend('匹配滤波', 'STFT替代', 'Location', 'southwest');
% end
% xlim([zoom_start, zoom_end]); ylim([-100, 5]); grid on;

% subplot(3,2,4);
% % 脉冲压缩效果时域对比
% t_comp = (0:199)*Ts*1e9; % 纳秒
% plot(t_comp, abs(Compressed_MF(1:200))/max(abs(Compressed_MF)), 'b-', 'LineWidth', 2);
% hold on;
% plot(t_comp, abs(Compressed_FrFT(1:200))/max(abs(Compressed_FrFT)), 'r-', 'LineWidth', 2);
% xlabel('时间 (ns)'); ylabel('归一化幅度');
% title('脉冲压缩效果对比');
% if FrFT_available
%     legend('匹配滤波', 'FrFT'); 
% else
%     legend('匹配滤波', 'STFT替代');
% end
% grid on;

% % 标记-3dB宽度
% half_power_level = 1/sqrt(2);
% plot([min(t_comp), max(t_comp)], [half_power_level, half_power_level], 'k--', 'LineWidth', 1);
% 
% subplot(3,2,5);
% % 动态范围对比
% methods = {'传统OTDR', '匹配滤波', 'FrFT/STFT'};
% dynamic_ranges = [dynamic_range_Trad, dynamic_range_MF, dynamic_range_FrFT];
% bar(dynamic_ranges, 'FaceColor', [0.7 0.7 0.9]);
% set(gca, 'XTickLabel', methods);
% ylabel('动态范围 (dB)');
% title('动态范围对比');
% grid on;
% 
% % 在柱状图上显示数值
% for i = 1:length(dynamic_ranges)
%     text(i, dynamic_ranges(i)+1, sprintf('%.1f dB', dynamic_ranges(i)), ...
%          'HorizontalAlignment', 'center', 'FontWeight', 'bold');
% end
% 
% subplot(3,2,6);
% % 分辨率对比
% resolution_methods = {'理论', '匹配滤波', 'FrFT/STFT'};
% resolutions = [theoretical_resolution, achieved_resolution_MF, achieved_resolution_FrFT];
% bar(resolutions, 'FaceColor', [0.9 0.7 0.7]);
% set(gca, 'XTickLabel', resolution_methods);
% ylabel('分辨率 (m)');
% title('分辨率对比');
% grid on;
% 
% % 在柱状图上显示数值
% for i = 1:length(resolutions)
%     text(i, resolutions(i)+2, sprintf('%.1f m', resolutions(i)), ...
%          'HorizontalAlignment', 'center', 'FontWeight', 'bold');
% end

% 图3: 最终OTDR轨迹与性能总结
figure('Position', [200, 200, 1200, 800], 'Name', '最终OTDR轨迹与性能总结');
plot(distance_axis/1000, Compressed_MF_dB, 'b-', 'LineWidth', 1.5);
hold on;
plot(distance_axis/1000, Traditional_dB, 'r-', 'LineWidth', 1.5);
% 添加LFM脉冲的动态范围标注
line([0, FiberLength/1000], [signal_MF, signal_MF], 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.5);
line([0, FiberLength/1000], [noise_floor_MF, noise_floor_MF], 'Color', 'g', 'LineStyle', ':', 'LineWidth', 1.5);
text(FiberLength/1000*0.7, signal_MF+5, sprintf('LFM初始散射功率: %.1f dB', signal_MF), 'Color', 'g', 'FontSize', 10);
text(FiberLength/1000*0.7, noise_floor_MF-8, sprintf('LFM噪声基底: %.1f dB', noise_floor_MF), 'Color', 'g', 'FontSize', 10);
text(FiberLength/1000*0.7, (signal_MF+noise_floor_MF)/2, sprintf('LFM动态范围: %.1f dB', dynamic_range_MF), 'Color', 'b', 'FontSize', 10, 'FontWeight', 'bold');

% 添加传统脉冲的动态范围标注
line([0, FiberLength/1000], [signal_TRA, signal_TRA], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);
line([0, FiberLength/1000], [noise_floor_Trad, noise_floor_Trad], 'Color', 'r', 'LineStyle', ':', 'LineWidth', 1.5);
text(FiberLength/1000*0.1, signal_TRA+5, sprintf('传统初始散射功率: %.1f dB', signal_TRA), 'Color', 'b', 'FontSize', 10);
text(FiberLength/1000*0.1, noise_floor_Trad-8, sprintf('传统噪声基底: %.1f dB', noise_floor_Trad), 'Color', 'b', 'FontSize', 10);
text(FiberLength/1000*0.1, (signal_TRA+noise_floor_Trad)/2, sprintf('传统动态范围: %.1f dB', dynamic_range_Trad), 'Color', 'b', 'FontSize', 10, 'FontWeight', 'bold');

xlabel('光纤距离 (km)'); ylabel('强度 (dB)');
title('200km LFM-OTDR vs 传统脉冲OTDR 轨迹对比');
for i = 1: length(time_idex_all)
n = time_idex_all(i); 
% 2. 提取该点的 X 和 Y 坐标
% 注意：必须保持和 plot 函数中完全一致的变换
x_val = distance_axis(n) / 1000 ;      % X轴坐标
y_val = Compressed_MF_dB(n);    % Y轴坐标

% 3. 锁定当前图层，防止覆盖
hold on; 

% 4. (可选) 在该点画一个红色的圆圈标记，使其更明显
plot(x_val, y_val, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);

% 5. 添加箭头和数值文本
% str 是显示的字符串：\leftarrow 生成左箭头，num2str 将数值转为字符串
str = [' \leftarrow ' num2str(y_val)]; 
text(x_val, y_val, str, 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'red');

% 6. 解除锁定
hold off;
end
for i = 1: length(time_idex_all)
n = time_idex_all(i); 
% 2. 提取该点的 X 和 Y 坐标
% 注意：必须保持和 plot 函数中完全一致的变换
x_val = distance_axis(n) / 1000 ;      % X轴坐标
y_val = Traditional_dB(n);    % Y轴坐标

% 3. 锁定当前图层，防止覆盖
hold on; 

% 4. (可选) 在该点画一个红色的圆圈标记，使其更明显
plot(x_val, y_val, 'go', 'MarkerSize', 8, 'LineWidth', 1.5);

% 5. 添加箭头和数值文本
% str 是显示的字符串：\leftarrow 生成左箭头，num2str 将数值转为字符串
str = [' \leftarrow ' num2str(y_val)]; 
text(x_val, y_val, str, 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'red');

% 6. 解除锁定
hold off;
end
xlim([0, (FiberLength-3000)/1000]); ylim([-150, -40]); 
grid on;

% 添加图例
legend('LFM-OTDR', '传统OTDR', 'LFM初始散射功率', 'LFM噪声基底', '传统初始散射功率', '传统噪声基底', ...
       'Location', 'best', 'NumColumns', 2);

% % 标记已知事件
% colors = ['g', 'm', 'c', 'y'];
% for i = 1:size(Events, 1)
%     event_km = Events(i, 1)/1000;
%     plot([event_km, event_km], [-100, 5], '--', 'Color', colors(i), 'LineWidth', 1);
%     text(event_km, -80 + mod(i,4)*8, sprintf('E%d@%.0fkm', i, event_km), ...
%          'FontSize', 9, 'Color', colors(i), 'BackgroundColor', 'white');
% end
% legend('LFM-OTDR(匹配滤波)', '传统OTDR', 'Location', 'southeast');
% 
% subplot(2,1,2);
% % 性能总结表格
% performance_metrics = {
%     '参数', '传统OTDR', 'LFM-OTDR(匹配滤波)', 'LFM-OTDR(FrFT/STFT)';
%     '动态范围(dB)', sprintf('%.1f', dynamic_range_Trad), sprintf('%.1f', dynamic_range_MF), sprintf('%.1f', dynamic_range_FrFT);
%     '分辨率(m)', '75.0', sprintf('%.1f', achieved_resolution_MF), sprintf('%.1f', achieved_resolution_FrFT);
%     '探测距离(km)', '200', '200', '200';
%     '事件检测数', sprintf('%d/%d', sum(Events(:,2) > noise_floor_Trad+3), size(Events,1)), ...
%                   sprintf('%d/%d', sum(Events(:,2) > noise_floor_MF+3), size(Events,1)), ...
%                   sprintf('%d/%d', sum(Events(:,2) > noise_floor_FrFT+3), size(Events,1));
%     '信噪比改善(dB)', '0.0', sprintf('+%.1f', dynamic_range_MF - dynamic_range_Trad), ...
%                      sprintf('+%.1f', dynamic_range_FrFT - dynamic_range_Trad);
% };
% 
% % 创建表格
% axis off;
% table_position = [0.1, 0.1, 0.8, 0.8];
% uitable('Data', performance_metrics, 'Position', table_position, ...
%         'ColumnWidth', {150, 120, 150, 150}, 'FontSize', 12, ...
%         'RowName', [], 'ColumnName', []);
% 
% title_handle = title('LFM-OTDR系统性能总结', 'FontSize', 14, 'FontWeight', 'bold');
% set(title_handle, 'Position', [0.5, 0.95, 0]);

fprintf('\n=== 仿真完成 ===\n');
fprintf('关键结论:\n');
fprintf('  - 理论分辨率: %.1f m\n', theoretical_resolution);

fprintf('  - LFM_OTDR: 动态范围=%.1f dB，初始散射功率=%.1f dB,噪声基底=%.1f dB，压缩后分辨率=%.1f m， 压缩前分辨率=%.1f m\n', dynamic_range_MF,signal_MF,noise_floor_MF, achieved_resolution_MF, achieved_resolution_MF1);
% if FrFT_available
%     fprintf('  - FrFT方法: 动态范围=%.1f dB, 分辨率=%.1f m\n', dynamic_range_FrFT, achieved_resolution_FrFT);
% else
%     fprintf('  - STFT替代方法: 动态范围=%.1f dB, 分辨率=%.1f m\n', dynamic_range_FrFT, achieved_resolution_FrFT);
% end
fprintf('  - 传统OTDR: 动态范围=%.1f dB,初始散射功率=%.1f dB,噪声基底=%.1f dB,分辨率:=%.1f m\n', dynamic_range_Trad,signal_TRA,noise_floor_Trad, achieved_resolution_tra);
% fprintf('  - 调制非线性误差: 最大约%.1f%%\n', max(abs(nonlinear_error)));

% %% ====================== 改进的STFT压缩函数 ======================
% function compressed_signal = improved_stft_compression(received_signal, reference_signal, Fs, mu)
%     % 改进的短时傅里叶变换压缩方法 (FrFT替代方案)
%     
%     N = length(received_signal);
%     M = length(reference_signal);
%     
%     % 使用STFT分析信号
%     window_length = min(512, M);
%     noverlap = round(window_length * 0.75);
%     nfft = max(1024, window_length);
%     
%     % 对接收信号进行STFT
%     [S_received, F_received, T_received] = spectrogram(received_signal, ...
%         window_length, noverlap, nfft, Fs, 'yaxis');
%     
%     % 对参考信号进行STFT
%     reference_padded = [reference_signal; zeros(N-M, 1)];
%     [S_reference, F_reference, T_reference] = spectrogram(reference_padded, ...
%         window_length, noverlap, nfft, Fs, 'yaxis');
%     
%     % 在时频域进行相关 (简化版本)
%     S_correlated = S_received .* conj(S_reference);
%     
%     % 逆STFT重建信号
%     compressed_signal = istft(S_correlated, Fs, 'Window', hamming(window_length), ...
%         'OverlapLength', noverlap, 'FFTLength', nfft);
%     
%     % 确保输出长度匹配
%     if length(compressed_signal) > N
%         compressed_signal = compressed_signal(1:N);
%     elseif length(compressed_signal) < N
%         compressed_signal = [compressed_signal; zeros(N-length(compressed_signal), 1)];
%     end
% end
