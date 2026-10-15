"""卡爾曼濾波＋RTS 平滑：由編碼器角度與馬達電流估計速度、加速度與「模型解釋不了的力矩」d(t)。

運動學（只用編碼器）  kin_rts
  狀態 [q, v, a]，白色急跳度（white jerk）模型，急跳度功率譜密度 sig_j²；量測 q（雜訊 sig_q）
  sig_q、sig_j 以新息對數概似最大化調整（tune_kin）。實測 sig_q ≈ 2e-5 rad（約 2 個編碼器 LSB），
  比純量化誤差（LSB/√12 ≈ 3.4e-6）大；若把 sig_q 固定在量化誤差，濾波器會把角度雜訊當成運動，
  估出的加速度雜訊很大。

加入力矩（擴增狀態）  kf_rts
  狀態 x = [q, v, a, d]，d 為隨機漫步（功率譜密度 sig_d²），即名義模型以外的力矩（殘差）
  量測 z1 = q；z2 = τ − known = hv·v + J·a + d（雜訊 sig_tau，只在 use_tau 為真的樣本使用）
       名義三參數模型：hv = B（依方向），known = Tc·sign(指令速度)
       擴充模型時，known 再加入已鑑別的漣波、溫度、行程、加減速階段等項
  起步、停止的暫態與極低速不在模型適用範圍，由呼叫端把 use_tau 設為假；此時 d 只隨機漫步。
  sig_q、sig_j 沿用運動學的調整值，只調 sig_tau、sig_d（tune_noise）。
"""
import numpy as np

DT = 1 / 125.0


def _model_mats(dt, sig_j, sig_d):
    A = np.array([[1, dt, dt ** 2 / 2, 0], [0, 1, dt, 0], [0, 0, 1, 0], [0, 0, 0, 1.0]])
    Qj = sig_j ** 2 * np.array([[dt ** 5 / 20, dt ** 4 / 8, dt ** 3 / 6],
                                [dt ** 4 / 8, dt ** 3 / 3, dt ** 2 / 2],
                                [dt ** 3 / 6, dt ** 2 / 2, dt]])
    Q = np.zeros((4, 4)); Q[:3, :3] = Qj; Q[3, 3] = sig_d ** 2 * dt
    return A, Q


def kin_rts(q, sig_q, sig_j, dt=DT, smooth=True):
    """只用編碼器：回傳 dict(xs = 平滑狀態 n×3 [q, v, a]（smooth 時）, loglik)。"""
    A4, Q4 = _model_mats(dt, sig_j, 0.0); A, Q = A4[:3, :3], Q4[:3, :3]
    r = sig_q ** 2; n = len(q)
    x = np.array([q[0], 0.0, 0.0]); P = np.diag([1e-8, 1e-2, 1.0])
    xp = np.empty((n, 3)); Pp = np.empty((n, 3, 3)); xf = np.empty((n, 3)); Pf = np.empty((n, 3, 3)); ll = 0.0
    for k in range(n):
        if k:
            x = A @ x; P = A @ P @ A.T + Q
        xp[k] = x; Pp[k] = P
        S = P[0, 0] + r; K = P[:, 0] / S; y = q[k] - x[0]
        x = x + K * y; P = P - np.outer(K, P[0, :])
        xf[k] = x; Pf[k] = P
        ll -= 0.5 * (np.log(2 * np.pi * S) + y * y / S)
    out = dict(loglik=ll)
    if smooth:
        xs = xf.copy()
        for k in range(n - 2, -1, -1):
            C = Pf[k] @ A.T @ np.linalg.inv(Pp[k + 1])
            xs[k] = xf[k] + C @ (xs[k + 1] - xp[k + 1])
        out["xs"] = xs
    return out


def tune_kin(q, init=(2e-5, 1.0), dt=DT):
    """以新息對數概似最大化調整 sig_q、sig_j（log 尺度，Nelder–Mead）。"""
    from scipy.optimize import minimize
    r = minimize(lambda th: -kin_rts(q, np.exp(th[0]), np.exp(th[1]), dt, smooth=False)["loglik"], np.log(init),
                 method="Nelder-Mead", options=dict(xatol=0.02, fatol=0.5, maxiter=200))
    return dict(sig_q=float(np.exp(r.x[0])), sig_j=float(np.exp(r.x[1])))


def kf_rts(q, z_tau, use_tau, hv, J, noise, dt=DT, smooth=True):
    """回傳 dict(xs = 平滑狀態 n×4, xf = 濾波狀態, nu_tau = 電流通道新息, S_tau = 新息變異, loglik)。
    noise = dict(sig_q, sig_tau, sig_j, sig_d)。"""
    n = len(q)
    A, Q = _model_mats(dt, noise["sig_j"], noise["sig_d"])
    rq, rt = noise["sig_q"] ** 2, noise["sig_tau"] ** 2
    x = np.array([q[0], 0.0, 0.0, 0.0]); P = np.diag([1e-8, 1e-2, 1.0, 100.0])
    xp = np.empty((n, 4)); Pp = np.empty((n, 4, 4)); xf = np.empty((n, 4)); Pf = np.empty((n, 4, 4))
    nu_t = np.full(n, np.nan); S_t = np.full(n, np.nan); ll = 0.0
    H1 = np.array([[1.0, 0, 0, 0]]); R1 = np.array([[rq]])
    for k in range(n):
        if k:
            x = A @ x; P = A @ P @ A.T + Q
        xp[k] = x; Pp[k] = P
        if use_tau[k]:
            H = np.array([[1.0, 0, 0, 0], [0, hv[k], J, 1.0]]); z = np.array([q[k], z_tau[k]]); R = np.diag([rq, rt])
        else:
            H = H1; z = np.array([q[k]]); R = R1
        y = z - H @ x
        S = H @ P @ H.T + R
        Si = np.linalg.inv(S)
        K = P @ H.T @ Si
        x = x + K @ y
        P = P - K @ S @ K.T; P = (P + P.T) / 2
        xf[k] = x; Pf[k] = P
        ll -= 0.5 * (len(y) * np.log(2 * np.pi) + np.log(np.linalg.det(S)) + y @ Si @ y)
        if use_tau[k]:
            nu_t[k] = y[1]; S_t[k] = S[1, 1]
    out = dict(xf=xf, nu_tau=nu_t, S_tau=S_t, loglik=ll)
    if smooth:
        xs = xf.copy()
        for k in range(n - 2, -1, -1):
            C = Pf[k] @ A.T @ np.linalg.inv(Pp[k + 1])
            xs[k] = xf[k] + C @ (xs[k + 1] - xp[k + 1])
        out["xs"] = xs
    return out


def tune_noise(q, z_tau, use_tau, hv, J, init, fixed=("sig_q", "sig_j"), dt=DT):
    """以新息對數概似最大化調整雜訊參數（只調 init 中不在 fixed 的項，log 尺度，Nelder–Mead）。"""
    from scipy.optimize import minimize
    keys = [k for k in init if k not in fixed]

    def nll(th):
        nz = dict(init); nz.update({k: float(np.exp(v)) for k, v in zip(keys, th)})
        return -kf_rts(q, z_tau, use_tau, hv, J, nz, dt, smooth=False)["loglik"]

    r = minimize(nll, np.log([init[k] for k in keys]), method="Nelder-Mead",
                 options=dict(xatol=0.02, fatol=0.5, maxiter=150))
    best = dict(init); best.update({k: float(np.exp(v)) for k, v in zip(keys, r.x)})
    return best, r
