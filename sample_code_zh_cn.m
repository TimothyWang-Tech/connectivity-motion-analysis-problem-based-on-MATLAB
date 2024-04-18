% 参数设置（无重力加速度有初速度情况）
Ni = 61; % 网格数
Nt = 154; % 时间步数
dksi = 1 / (Ni - 1); % 空间步长
dt = 0.01; % 时间步长
lambda = 1; % 常数
g = 0; % 重力加速度
v0 = 0.5 * sqrt(3); % 初始速度

% 初始化变量
x = zeros(Ni, Nt); % x坐标
y = zeros(Ni, Nt); % y坐标
T = zeros(Ni, Nt); 
p = zeros(1, Nt); % 界面位置

% j = 0 时的初始化
p(1) = 0.5 * (Ni + 1);
for i = 1 : p(1) 
    x(i, 1) = 0;
    y(i, 1) = -0.5 + (i - 1) * dksi;
    T(i, 1) = 0.5 * lambda * g * (i - 1) * dksi;
end

% j = 1 时的初始化
p(2) = p(1) + v0 * dt / dksi;
for i = 1 : p(1) - 1
    x(i, 2) = 0;
    y(i, 2) = y(i, 1) - v0 * dt;
end
x(fix(p(1)), 2) = v0 * dt;
y(fix(p(1)), 2) = y(fix(p(1)), 1);
x(fix(p(1)) + 1, 2) = (p(1) - (fix(p(1)) + 1)) * dksi;

% j = 1 时的求解
for j = 2 : Nt
    % j 时刻的矩阵维数
    n = fix(p(j));
    % 创建稀疏矩阵 A 和向量 b
    A = spdiags(zeros(n, 3), [-1, 0, 1], n - 1, n - 1);
    b = zeros(n - 1, 1);
    % 填充 A 和 b 的元素
    for i = 3 : n
        youcexishu = -(x(i, j) - x(i - 1, j)) * (x(i, j) - x(i - 1, j) - x(i, j - 1) + x(i - 1, j - 1))- (y(i, j) - y(i - 1, j)) * (y(i, j) - y(i - 1, j) - y(i, j - 1) + y(i - 1, j - 1));
        b(i - 2) = youcexishu;
    end
    b(n - 1) = -(x(n + 1, j) - x(n, j)) * (dksi * (p(j) - p(j - 1)) - x(n, j) + x(n, j - 1)) - y(n, j) * (y(n, j) - y(n, j - 1));
    b = sparse(b);
    for i = 3 : 3
       A(i - 2, i - 2) = (dt * dt / (dksi * dksi)) * (2 * (x(i - 1, j) - x(i, j)) * (x(i, j) - x(i - 1, j)) + 2 * (y(i - 1, j) - y(i, j)) * (y(i, j) - y(i - 1, j)));
        A(i - 2, i - 1) = (dt * dt / (dksi * dksi)) * ((x(i + 1, j) - x(i, j)) * (x(i, j) - x(i - 1, j)) + (y(i + 1, j) - y(i, j)) * (y(i, j) - y(i - 1, j)));
        A(i - 2, i - 2) = A(i - 2, i - 2);
        A(i - 2, i - 1) = A(i - 2, i - 1);
    end
    for i = 4 : n
        A(i - 2, i - 3) = -(dt * dt / (dksi * dksi)) * ((x(i - 2, j) - x(i - 1, j)) * (x(i, j) - x(i - 1, j)) + (y(i - 2, j) - y(i - 1, j)) * (y(i, j) - y(i - 1, j)));
        A(i - 2, i - 2) = (dt * dt / (dksi * dksi)) * (2 * (x(i - 1, j) - x(i, j)) * (x(i, j) - x(i - 1, j)) + 2 * (y(i - 1, j) - y(i, j)) * (y(i, j) - y(i - 1, j)));
        A(i - 2, i - 1) = (dt * dt / (dksi * dksi)) * ((x(i + 1, j) - x(i, j)) * (x(i, j) - x(i - 1, j)) + (y(i + 1, j) - y(i, j)) * (y(i, j) - y(i - 1, j)));
        A(i - 2, i - 2) = A(i - 2, i - 2);
        A(i - 2, i - 1) = A(i - 2, i - 1);
        A(i - 2, i - 3) = A(i - 2, i - 3);
    end

    A(n - 1, n - 2) = (dt * dt / (dksi * dksi)) * ((x(n + 1, j) - x(n, j)) * (x(n - 1, j) - x(n, j)) + y(n, j) * (y(n - 1, j) - y(n, j)));
    A(n - 1, n - 1) = (dt * dt / (dksi * dksi)) * ((x(n + 1, j) - x(n, j)) * ((dksi * dksi / (1 - p(j) * dksi)) - (x(n + 1, j) - x(n, j))) + y(n, j) * (y(n + 1, j) - y(n, j)));
    k = sprank(A);
    vecsovle = A \ b;
    vecsolvefull = full(vecsovle);
    for i = 2 : n
        T(i, j) = double(vecsolvefull(i - 1));
    end
    T(1, j) = 0;
    p(j + 1) = ((dt * dt * T(n, j)) / (dksi * lambda * (1 - p(j) * dksi))) + 2 * p(j) - p(j - 1);
    for i = 2 : n - 1
        x(i, j + 1) = (dt * dt / (dksi * dksi)) * (T(i - 1, j) * (x(i - 1, j) - x(i, j)) + T(i, j) * (x(i + 1, j) - x(i, j))) + 2 * x(i, j) - x(i, j - 1);
        y(i, j + 1) = (dt * dt / (dksi * dksi)) * (T(i - 1, j) * (y(i - 1, j) - y(i, j)) + T(i, j) * (y(i + 1, j) - y(i, j))) + 2 * y(i, j) - y(i, j - 1);
    end

    x(1, j + 1) = ((dt * dt) / dksi) * T(2, j) * (x(2, j) - x(1, j)) + 2 * x(1, j) - x(1, j - 1);
    y(1, j + 1) = ((dt * dt) / dksi) * T(2, j) * (y(2, j) - y(1, j)) + 2 * y(1, j) - y(1, j - 1);
    
    if fix(p(j)) == fix(p(j + 1))
        x(n, j + 1) = (dt * dt / (dksi * dksi)) * (T(n - 1, j) * (x(n - 1, j) - x(n, j)) + T(n, j) * (x(n + 1, j) - x(n, j))) + 2 * x(n, j) - x(n, j - 1);
        y(n, j + 1) = (dt * dt / (dksi * dksi)) * (T(n - 1, j) * (y(n - 1, j) - y(n, j)) + T(n, j) * (y(n + 1, j) - y(n, j))) + 2 * y(n, j) - y(n, j - 1);
        x(n + 1, j + 1) = (p(j + 1) - n - 1) * dksi;
        y(n + 1, j + 1) = 0;
    end
    
    if fix(p(j)) + 1 == fix(p(j + 1))
        x(n, j + 1) = (dt * dt / (dksi * dksi)) * (T(n - 1, j) * (x(n - 1, j) - x(n, j)) + T(n, j) * (x(n + 1, j) - x(n, j))) + 2 * x(n, j) - x(n, j - 1);
        y(n, j + 1) = (dt * dt / (dksi * dksi)) * (T(n - 1, j) * (y(n - 1, j) - y(n, j)) + T(n, j) * (y(n + 1, j) - y(n, j))) + 2 * y(n, j) - y(n, j - 1);
        x(n + 1, j + 1) = (p(j + 1) - n - 1) * dksi;
        y(n + 1, j + 1) = 0;
    end
end

for j = 1 : Nt + 1
    disp(['j=', num2str(j)]);
    disp(['p[j]=', num2str(p(j))]);
    for i = 1 : p(j) + 1
        fprintf('%d ', x(i, j));
    end
    disp(' ');
    for i = 1 : p(j) + 1
        fprintf('%d ', y(i, j));
    end
    disp(' ');
end

F = struct('cdata', [], 'colormap', []);

for j = 1 : Nt
    x0 = -(1 - p(j) * dksi);
    y0 = 0;
    plinex1 = [0, x0];
    pliney1 = [0, y0];
    plinex2 = [0, x(fix(p(j)), j)];
    pliney2 = [0, y(fix(p(j)), j)];
    plot(10 * x(2:p(j) + 1, j), 10 * y(2:p(j) + 1, j), 'o', 'Color', 'w', 'MarkerSize', 10);
    hold on;
    plot(x(2:p(j), j), y(2:p(j), j), '-', 'Color', 'b', 'LineWidth', 2);
    plot(plinex1, pliney1, '-', 'Color', 'b', 'LineWidth', 2);
    plot(plinex2, pliney2, '-', 'Color', 'b', 'LineWidth', 2);
    axis on
    axis equal
    set(gcf, 'color', 'w')
    xlim([-0.5 0.5])
    ylim([-1.3 0])
    F(j) = getframe;
    hold off;
end
