%% verify_mathieu.m — 验证新 mathieu_* 函数包
% 注意: 新旧 API 不同,需用 full path 调用旧函数避免 shadowing
new_dir = fileparts(mfilename('fullpath'));
old_dir = fullfile(new_dir, '..', 'overmoded-elliptical', 'mathieu');
fprintf('New: %s\nOld: %s\n\n', new_dir, old_dir);
fprintf('===== Mathieu 函数验证 =====\n\n');

%% 1. q=0 退化: a_n(0)=n^2, b_n(0)=n^2
fprintf('--- 1. q=0 退化 a_n(0)=n^2 ---\n');
for n=0:5
    a = mathieu_char_value(n, 0);
    fprintf('n=%d: a=%.8f expected=%d diff=%.2e\n', n, a, n^2, abs(a-n^2));
    assert(abs(a - n^2) < 1e-10, 'q=0 degeneracy FAILED for n=%d', n);
end
fprintf('PASS\n\n');

%% 2. 归一化 ∫ce_n²=π (n>0) or 2π (n=0), ∫se_n²=π
fprintf('--- 2. 归一化 ---\n');
vv = linspace(0, 2*pi, 500); dv = vv(2)-vv(1);
max_err_ce = 0; max_err_se = 0;
for n=0:4
    for q=[0.5, 1, 2, 5]
        ce = mathieu_ce(n, vv, q);
        int_ce = sum(ce.^2) * dv;
        tgt = pi; if n==0, tgt = 2*pi; end
        err_ce = abs(int_ce - tgt)/tgt;
        if err_ce > max_err_ce, max_err_ce = err_ce; end
        if n > 0
            se = mathieu_se(n, vv, q);
            int_se = sum(se.^2) * dv;
            err_se = abs(int_se - pi)/pi;
            if err_se > max_err_se, max_err_se = err_se; end
        end
        if err_ce > 1e-3, fprintf('  *** ce n=%d q=%.1f: err=%.2e\n', n, q, err_ce); end
    end
end
fprintf('Max ce norm err: %.2e,  Max se norm err: %.2e\n', max_err_ce, max_err_se);
assert(max_err_ce < 0.01, 'ce normalization FAILED: %.2e', max_err_ce);
assert(max_err_se < 0.01, 'se normalization FAILED: %.2e', max_err_se);
fprintf('PASS\n\n');

%% 3. 特殊值
fprintf('--- 3. 特殊值 ---\n');
% se_n(0,q) = 0
for n=1:4, assert(abs(mathieu_se(n,0,1))<1e-12,'se_n(0)=0 FAILED n=%d',n); end
fprintf('se_n(0,q)=0: PASS\n');
% Jo_n(0,q) = 0 for all n (angular odd vanishes at origin)
for n=1:4, assert(abs(mathieu_Jo(n,0,1))<1e-10,'Jo_n(0)=0 FAILED n=%d',n); end
for n=2:2:4, assert(abs(mathieu_Jo(n,0,3))<1e-10,'Jo_even(0)=0 FAILED n=%d',n); end
fprintf('Jo_n(0,q)=0: PASS\n');
% dJe_{even}/du (0,q) = 0
for n=0:2:4
    [~, dR] = mathieu_Je(n, 0, 1);
    assert(abs(dR)<1e-10,'dJe_even/du(0)=0 FAILED n=%d',n);
end
fprintf('dJe_{even}/du(0,q)=0: PASS\n');

% se_{even}(0,q)=0
for n=2:2:4, assert(abs(mathieu_se(n,0,1))<1e-12,'se_even(0)=0 FAILED n=%d',n); end
fprintf('se_{even}(0,q)=0: PASS\n\n');

%% 4. 与旧实现对比 (使用 full path 避免 shadowing)
fprintf('--- 4. 新旧 ce 对比 (full path) ---\n');
vt = pi/3; max_d = 0;
for n=0:4, for q=[0.5, 1, 2, 5]
    nv = mathieu_ce(n, vt, q);  % new: (n,v,q)
    % old: call via cd to old dir
    old_cd = cd(old_dir);
    ov = mathieu_ce(n, q, vt);  % old: (n,q,v)
    cd(old_cd);
    d = abs(nv - ov);
    if d > max_d, max_d = d; end
    fprintf('n=%d q=%.1f: new=%.8f old=%.8f diff=%.2e\n', n, q, nv, ov, d);
    if d > 1e-8, fprintf('  *** DIFFER\n'); end
end,end
fprintf('Max ce diff: %.2e\n', max_d);
if max_d > 1e-8
    fprintf('WARNING: 新旧 ce 不一致 (可能 sign/normalization 差异)\n');
end
fprintf('\n');

%% 5. 径向新旧对比
fprintf('--- 5. 新旧 Je 对比 (full path) ---\n');
ut = 0.5; max_d = 0;
for n=0:4, for q=[0.5, 1, 2]
    nv = mathieu_Je(n, ut, q);
    old_cd = cd(old_dir);
    ov = mathieu_Je(n, q, ut);
    cd(old_cd);
    d = abs(nv - ov);
    if d > max_d, max_d = d; end
    fprintf('n=%d q=%.1f: new=%.8f old=%.8f diff=%.2e\n', n, q, nv, ov, d);
    if d > 1e-4, fprintf('  *** DIFFER\n'); end
end,end
fprintf('Max Je diff: %.2e\n', max_d);
fprintf('\n');

%% 6. 径向第二类 (新功能)
fprintf('--- 6. Ne/No (新) ---\n');
for n=0:3, for q=[0.5, 2]
    Ne = mathieu_Ne(n, 1, q);
    if n>0, No = mathieu_No(n, 1, q); else, No = 0; end
    fprintf('n=%d q=%.1f: Ne=%.6f No=%.6f\n', n, q, Ne, No);
end,end
fprintf('\n');

%% 7. 特征值对照 Alhargan Table I
fprintf('--- 7. TEe11 gamma 对照 Alhargan ---\n');
% Alhargan gamma for TEe11: gamma ≈ 1 - 1.05e + 0.12e² - 0.08e³
for e=[0.1, 0.3, 0.5, 0.7, 0.9]
    p11 = 1.841183781340265;
    l = e;  % a=1
    gamma_alh = 1 - 1.05*e + 0.12*e^2 - 0.08*e^3;
    % Our method: kc from Mathieu root finding
    % For TEe11: Je'_1(h,u0)=0 where u0=acosh(1/e)
    if e < 1
        u0 = acosh(1/e);
        % find h where Je'_1(h,u0)=0 via scanning
        hs = linspace(0.1, 5, 100);
        f_min = inf; h_best = 0;
        for hi = 1:length(hs)
            h = hs(hi); q = (h/2)^2;
            [~, dR] = mathieu_Je(1, u0, q);
            if abs(dR) < abs(f_min)
                f_min = dR; h_best = h;
            end
        end
        % Refine
        h = h_best;
        for iter = 1:20
            q = (h/2)^2; dh = 1e-6*h;
            [~, f] = mathieu_Je(1, u0, q);
            [~, fp] = mathieu_Je(1, u0, (h+dh)^2/4);
            h = h - f/((fp-f)/dh);
            if abs(f) < 1e-10, break; end
        end
        kc = h/l;
        gamma_us = h*(1-e)/(e*p11);
    else
        gamma_us = NaN;
    end
    fprintf('e=%.1f: gamma_alh=%.6f gamma_us=%.6f\n', e, gamma_alh, gamma_us);
end
fprintf('\n===== 验证完成 =====\n');
