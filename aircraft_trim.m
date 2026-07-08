% 固定翼飞机配平：设定迎角定直平飞
clc; 
clear; 
close all;

%% 1. 气动数据
% 所有角度单位为度，气动系数基于度进行插值
aero = struct();

% 纵向基本曲线（迎角范围：-5° ~ 45°）
aero.aLon.alpha = (-5:2:45)';
aero.aLon.CL    = 0.2 + 0.1 * aero.aLon.alpha;          % 升力线斜率 ~0.1/deg
aero.aLon.CD    = 0.02 + 0.002 * (aero.aLon.alpha/10).^2; % 阻力曲线
aero.aLon.Cm    = -0.05 - 0.02 * aero.aLon.alpha;       % 俯仰力矩（静稳定）

% 升降舵舵效曲线（假设与迎角无关，常值导数）
aero.aDe.alpha  = aero.aLon.alpha;
aero.aDe.CLDe   = 0.01 * ones(size(aero.aDe.alpha));   % dCL/dδe
aero.aDe.CDDe   = 0.002 * ones(size(aero.aDe.alpha));  % dCD/dδe
aero.aDe.CmDe   = -0.03 * ones(size(aero.aDe.alpha));  % dCm/dδe（低头为负）

% 侧向基本曲线（侧滑角范围：-20° ~ 20°）
aero.aLat.beta  = (-20:5:20)';
aero.aLat.CY    = -0.01 * aero.aLat.beta;      % 侧力系数
aero.aLat.Cl    = -0.002 * aero.aLat.beta;     % 滚转力矩系数
aero.aLat.Cn    = 0.001 * aero.aLat.beta;      % 偏航力矩系数

% 副翼/方向舵舵效（常数）
aero.aAil.dCl_da = -0.01;   % 滚转力矩/副翼偏度
aero.aRud.dCn_dr = -0.005;  % 偏航力矩/方向舵偏度
aero.aRud.dCY_dr = 0.002;   % 侧力/方向舵偏度

fprintf('气动示例数据已构建，迎角范围 %.0f° ~ %.0f°\n', ...
    aero.aLon.alpha(1), aero.aLon.alpha(end));

%% 2. 设置飞机参数
fprintf('\n=== 设置飞机参数 ===\n');
aircraft.mass = 7000;   % kg
aircraft.g    = 9.81;
aircraft.S    = 25.0;       % 机翼面积 (m²)
aircraft.ba   = 3.0;      % 平均气动弦长 (m)
aircraft.aL   = 8.0;       % 展长 (m)
aircraft.W    = aircraft.mass * aircraft.g;

fprintf('  质量: %.1f kg\n', aircraft.mass);
fprintf('  重量: %.1f N\n', aircraft.W);
fprintf('  机翼面积: %.1f m²\n', aircraft.S);

%% 3. 设置配平条件和牛顿法参数
fprintf('\n=== 设置配平条件 ===\n');
alpha_deg = 5;                     % 目标迎角 (度)
alpha     = deg2rad(alpha_deg);     % 用于三角函数计算的弧度值
rho       = 1.225;                  % 海平面密度 (kg/m³)

% 其他飞行状态（定直平飞）
beta = 0;        % 无侧滑 (rad)
wx = 0; wy = 0; wz = 0;            % 无角速度 (rad/s)
deltaa  = 0;        % 副翼中立 (deg)
deltar  = 0;        % 方向舵中立 (deg)

% 控制面范围限制 (单位：度)
deltae_min = -30;
deltae_max = 30;
T_min      = 0;
T_max      = 200000;   % N

% 牛顿法参数
max_iter       = 100;
tol            = 1e-6;
relaxation     = 1.0;
damping_factor = 0.8;

% 有限差分步长
h_V      = 0.1;    % 速度扰动 (m/s)
h_deltae = 0.1;    % 升降舵扰动 (度)
h_T      = 10;     % 推力扰动 (N)

%% 4. 牛顿迭代法配平（三个方程，三个未知数：V, deltae, T）
fprintf('\n=== 开始牛顿迭代法配平 ===\n');
x0 = [100; 5; 10000];   % 初始值[V (m/s); deltae (deg); T (N)]
x  = x0;

% 历史记录预分配
history = struct();
history.V         = zeros(max_iter,1);
history.deltae    = zeros(max_iter,1);
history.T         = zeros(max_iter,1);
history.residual1 = zeros(max_iter,1);   % 升力残差
history.residual2 = zeros(max_iter,1);   % 力矩残差
history.residual3 = zeros(max_iter,1);   % 阻力残差
history.normR     = zeros(max_iter,1);
history.L_aero    = zeros(max_iter,1);
history.D         = zeros(max_iter,1);
history.Mz        = zeros(max_iter,1);

fprintf('迭代 |  速度(m/s) | 升降舵(°) |   推力(kN)   | 升力残差 | 力矩残差 | 阻力残差 | 残差范数\n');
fprintf('-------------------------------------------------------------------------------------------\n');

converged = false;
for iter = 1:max_iter
    V      = x(1);
    deltae = x(2);
    T      = x(3);

    % 计算当前残差
    [R, L_aero, D_val, Mz_val] = compute_residuals_full(V, deltae, T, alpha_deg, ...
        alpha, aircraft, rho, beta, wx, wy, wz, deltaa, deltar, aero);

    % 记录历史
    history.V(iter)         = V;
    history.deltae(iter)    = deltae;
    history.T(iter)         = T;
    history.residual1(iter) = R(1);
    history.residual2(iter) = R(2);
    history.residual3(iter) = R(3);
    history.normR(iter)     = norm(R);
    history.L_aero(iter)    = L_aero;
    history.D(iter)         = D_val;
    history.Mz(iter)        = Mz_val;

    % 显示迭代结果
    fprintf('%3d | %10.3f | %10.3f | %12.3f | %9.3e | %9.3e | %9.3e | %9.3e\n', ...
        iter, V, deltae, T/1000, R(1), R(2), R(3), norm(R));

    % 检查收敛
    if norm(R) < tol
        fprintf('\n 收敛成功！\n');
        converged = true;
        break;
    end

    % 计算雅可比矩阵（中心差分，3x3）
    J = zeros(3,3);
    [Rp, ~, ~, ~] = compute_residuals_full(V+h_V, deltae, T, alpha_deg, alpha, ...
        aircraft, rho, beta, wx, wy, wz, deltaa, deltar, aero);
    [Rm, ~, ~, ~] = compute_residuals_full(V-h_V, deltae, T, alpha_deg, alpha, ...
        aircraft, rho, beta, wx, wy, wz, deltaa, deltar, aero);
    J(:,1) = (Rp - Rm) / (2*h_V);

    [Rp, ~, ~, ~] = compute_residuals_full(V, deltae+h_deltae, T, alpha_deg, alpha, ...
        aircraft, rho, beta, wx, wy, wz, deltaa, deltar, aero);
    [Rm, ~, ~, ~] = compute_residuals_full(V, deltae-h_deltae, T, alpha_deg, alpha, ...
        aircraft, rho, beta, wx, wy, wz, deltaa, deltar, aero);
    J(:,2) = (Rp - Rm) / (2*h_deltae);

    [Rp, ~, ~, ~] = compute_residuals_full(V, deltae, T+h_T, alpha_deg, alpha, ...
        aircraft, rho, beta, wx, wy, wz, deltaa, deltar, aero);
    [Rm, ~, ~, ~] = compute_residuals_full(V, deltae, T-h_T, alpha_deg, alpha, ...
        aircraft, rho, beta, wx, wy, wz, deltaa, deltar, aero);
    J(:,3) = (Rp - Rm) / (2*h_T);

    % 牛顿步
    if rcond(J) > eps
        dx = J \ (-R);
    else
        dx = pinv(J) * (-R);
        fprintf('  警告：雅可比矩阵接近奇异，使用伪逆\n');
    end

    dx = dx * relaxation;
    history.step_size = zeros(max_iter,1);
    history.step_size(iter) = norm(dx);

    % 状态更新与边界限制
    x_new = x + dx;
    x_new(1) = max(30, min(300, x_new(1)));             % 速度限制 30~300 m/s
    x_new(2) = max(deltae_min, min(deltae_max, x_new(2))); % 升降舵限制 (度)
    x_new(3) = max(T_min, min(T_max, x_new(3)));        % 推力限制

    if x_new(1) == 30 || x_new(1) == 300
        fprintf('  警告：速度达到边界\n');
    end
    if x_new(2) == deltae_min || x_new(2) == deltae_max
        fprintf('  警告：升降舵达到边界\n');
    end
    if x_new(3) == T_min || x_new(3) == T_max
        fprintf('  警告：推力达到边界\n');
    end

    % 线搜索（如残差增大）
    [R_new, ~, ~, ~] = compute_residuals_full(x_new(1), x_new(2), x_new(3), alpha_deg, alpha, ...
        aircraft, rho, beta, wx, wy, wz, deltaa, deltar, aero);
    if norm(R_new) > norm(R) && iter > 1
        for ls = 1:6
            dx_scaled = dx * 0.5^ls;
            x_trial = x + dx_scaled;
            x_trial(1) = max(30, min(300, x_trial(1)));
            x_trial(2) = max(deltae_min, min(deltae_max, x_trial(2)));
            x_trial(3) = max(T_min, min(T_max, x_trial(3)));
            [R_trial, ~, ~, ~] = compute_residuals_full(x_trial(1), x_trial(2), x_trial(3), alpha_deg, alpha, ...
                aircraft, rho, beta, wx, wy, wz, deltaa, deltar, aero);
            if norm(R_trial) < norm(R)
                x_new = x_trial;
                fprintf('  线搜索：步长缩减为 %.4f\n', 0.5^ls);
                break;
            end
        end
    end

    x = x_new;

    % 阻尼因子自适应调整
    if iter > 1
        if iter > 2 && history.normR(iter) > 0.9*history.normR(iter-1)
            relaxation = max(relaxation * damping_factor, 0.1);
            fprintf('  阻尼因子调整为: %.2f\n', relaxation);
        elseif history.normR(iter) < 0.5*history.normR(iter-1)
            relaxation = min(relaxation / damping_factor, 1.0);
        end
    end
end

if ~converged
    fprintf('\n 达到最大迭代次数，未收敛\n');
    iter_actual = max_iter;
else
    iter_actual = iter;
end

% 裁剪历史记录
flds = fieldnames(history);
for i = 1:numel(flds)
    history.(flds{i}) = history.(flds{i})(1:iter_actual);
end

%% 5. 输出配平结果
fprintf('\n=== 配平结果 ===\n');
fprintf('目标迎角: %.1f°\n', alpha_deg);
fprintf('配平速度: %.6f m/s\n', x(1));
fprintf('配平升降舵偏度: %.6f°\n', x(2));
fprintf('所需推力: %.6f kN\n', x(3)/1000);
fprintf('迭代次数: %d\n', iter_actual);
fprintf('最终残差范数: %.3e\n', history.normR(end));

% 检查受力平衡
fprintf('\n力与力矩平衡验证:\n');
[~, L, D, Mz] = compute_residuals_full(x(1), x(2), x(3), alpha_deg, alpha, ...
    aircraft, rho, beta, wx, wy, wz, deltaa, deltar, aero);
fprintf('气动升力: %.2f N\n', L);
fprintf('重力：%.2f N \n',aircraft.W);
fprintf('推力: %.2f N\n', x(3));
fprintf('阻力：%.2 N\n',D);
fprintf('俯仰力矩: %.2f N·m\n', Mz);

% 简单的收敛曲线图
figure;
subplot(2,1,1);
semilogy(1:iter_actual, history.normR, 'b-o');
xlabel('迭代次数'); ylabel('残差范数'); grid on;
title('配平迭代收敛历程');

subplot(2,1,2);
yyaxis left;
plot(1:iter_actual, history.V, 'b-', 'LineWidth',1.5); ylabel('速度 (m/s)');
yyaxis right;
plot(1:iter_actual, history.T/1000, 'r-', 'LineWidth',1.5); ylabel('推力 (kN)');
xlabel('迭代次数'); grid on;
legend('速度','推力','Location','best');
title('状态量收敛历史');


% =========================================================================
function [residual, L_aero, D, Mz] = compute_residuals_full(V, deltae_deg, T, alpha_deg, alpha, aircraft, rho, ...
    beta, wx, wy, wz, deltaa_deg, deltar_deg, aero)
    % 计算配平残差（升力、俯仰力矩、阻力方向）
    % 输入：V (m/s), deltae_deg (度), T (N), alpha_deg (度, 用于气动插值),
    %       alpha (rad, 用于三角函数), 其余参数结构体
    % 输出：residual (3x1), L_aero (N), D (N), Mz (N·m)

    % 计算体轴系速度分量（Z轴向下为正）
    Vx = V * cos(alpha) * cos(beta);
    Vy = V * sin(beta);
    Vz = V * sin(alpha) * cos(beta);   % 负号表示来流在Z方向分量为负（向下为正）

    % 调用气动模型（输出包含体轴力和风轴力/力矩）
    [~, ~, ~, ~, ~, ~, L_aero, D, Mz] = aerodynamicloads(Vx, Vy, Vz, wx, wy, wz, ...
        deltae_deg, deltar_deg, deltaa_deg, aero, aircraft, rho);

    % 总升力（含推力分量）
    L_total = L_aero + T * sin(alpha);

    % 三个残差方程（无量纲化）
    residual = zeros(3,1);
    residual(1) = (L_total - aircraft.W) / aircraft.W;          % 升力平衡
    residual(2) = Mz / (aircraft.W * aircraft.ba);               % 俯仰力矩平衡
    residual(3) = (T * cos(alpha) - D) / aircraft.W;             % 阻力平衡
end

% -------------------------------------------------------------------------
function [X_body, Y_body, Z_body, Mx, My, Mz, L, D, M_pitch] = aerodynamicloads(...
    Vx, Vy, Vz, wx, wy, wz, delta_e_deg, delta_r_deg, delta_a_deg, aero, veh, rho)
    % 气动力/力矩计算（插值法，所有角度输入为度）
    % 输入：
    %   Vx,Vy,Vz - 体轴系速度分量 (m/s)，Z向下为正
    %   wx,wy,wz - 角速度 (rad/s)（本脚本未使用，预留）
    %   delta_e_deg, delta_r_deg, delta_a_deg - 舵面偏度 (度)
    % 输出：
    %   X_body, Y_body, Z_body - 体轴系力 (N)，向前/向右/向下
    %   Mx, My, Mz - 体轴系力矩 (N·m)
    %   L, D - 风轴系升力、阻力 (N)
    %   M_pitch - 俯仰力矩 (N·m)，即My

    % 计算空速和气流角
    V = sqrt(Vx^2 + Vy^2 + Vz^2);
    if V < 1e-6
        X_body=0; Y_body=0; Z_body=0; Mx=0; My=0; Mz=0;
        L=0; D=0; M_pitch=0;
        return;
    end
    alpha = atan2(Vz, Vx);          % 体轴系迎角 (rad)，Vz向下为正时，α>0 表示机头上仰
    beta  = asin(Vy / V);           % 侧滑角 (rad)

    % 转换为度用于插值
    alpha_deg = rad2deg(alpha);
    beta_deg  = rad2deg(beta);

    qdyn = 0.5 * rho * V^2;

    % --- 纵向系数插值 ---
    CL0 = interp1(aero.aLon.alpha, aero.aLon.CL, alpha_deg, 'linear', 'extrap');
    CD0 = interp1(aero.aLon.alpha, aero.aLon.CD, alpha_deg, 'linear', 'extrap');
    Cm0 = interp1(aero.aLon.alpha, aero.aLon.Cm, alpha_deg, 'linear', 'extrap');

    % 升降舵贡献（舵效导数乘偏度）
    CLde = interp1(aero.aDe.alpha, aero.aDe.CLDe, alpha_deg, 'linear', 'extrap') * delta_e_deg;
    CDde = interp1(aero.aDe.alpha, aero.aDe.CDDe, alpha_deg, 'linear', 'extrap') * delta_e_deg;
    Cmde = interp1(aero.aDe.alpha, aero.aDe.CmDe, alpha_deg, 'linear', 'extrap') * delta_e_deg;

    CL = CL0 + CLde;
    CD = CD0 + CDde;
    Cm = Cm0 + Cmde;

    % 风轴系升力和阻力
    L = qdyn * veh.S * CL;
    D = qdyn * veh.S * CD;

    % --- 侧向系数插值 ---
    CY0 = interp1(aero.aLat.beta, aero.aLat.CY, beta_deg, 'linear', 'extrap');
    CYdr = aero.aRud.dCY_dr * delta_r_deg;
    CY   = CY0 + CYdr;

    Cl0  = interp1(aero.aLat.beta, aero.aLat.Cl, beta_deg, 'linear', 'extrap');
    Clda = aero.aAil.dCl_da * delta_a_deg;
    Cl   = Cl0 + Clda;

    Cn0  = interp1(aero.aLat.beta, aero.aLat.Cn, beta_deg, 'linear', 'extrap');
    Cndr = aero.aRud.dCn_dr * delta_r_deg;
    Cn   = Cn0 + Cndr;

    % 侧力与力矩
    Y_body  = qdyn * veh.S * CY;
    Mx_roll = qdyn * veh.S * veh.aL * Cl;   % 滚转力矩 (绕X轴)
    Mz_yaw  = qdyn * veh.S * veh.aL * Cn;   % 偏航力矩 (绕Z轴)

    % --- 转换风轴力到体轴系 ---
    % 体轴系：X向前, Z向下
    X_body = -D * cos(alpha) + L * sin(alpha);
    Z_body = -D * sin(alpha) - L * cos(alpha);

    % 俯仰力矩（绕Y轴）
    My = qdyn * veh.S * veh.ba * Cm;
    Mx = Mx_roll;
    Mz = Mz_yaw;

    % 额外输出（风轴）
    M_pitch = My;   % 俯仰力矩
end