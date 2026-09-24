"""Smooth Huber regression on a fixed optical basis, with optional base weights.

QR/SVD preconditioning preserves the fitted function space. The intercept is
fitted explicitly. No normal equations of the original polynomial are formed.
"""
import numpy as np
from scipy.linalg import qr, svd, lstsq
from scipy.optimize import minimize


def loss_gradient(beta, q, target, delta, weights=None):
    residual = q @ beta - target
    norm = np.hypot(1., residual / delta)
    # Algebraically delta**2*(norm-1), avoiding cancellation near zero.
    weights = np.ones(len(target)) if weights is None else weights
    loss = np.sum(weights * residual**2 / (norm + 1.)) / weights.sum()
    gradient = q.T @ (weights * residual / norm) / weights.sum()
    return float(loss), gradient


class FixedBasis:
    def __init__(self, x, rcond=1e-12, weights=None):
        x = np.asarray(x, dtype=float)
        if x.ndim != 2 or len(x) <= x.shape[1] or not np.isfinite(x).all():
            raise ValueError('Require a finite, overdetermined design')
        if not np.array_equal(x[:, 0], np.ones(len(x))):
            raise ValueError('First column must be the explicit constant')
        w = np.ones(len(x)) if weights is None else np.asarray(weights, dtype=float)
        if w.shape != (len(x),) or not np.isfinite(w).all() or np.any(w <= 0):
            raise ValueError('Base weights must be finite, positive and match the rows')
        self.weights = w / w.mean()
        self.mean = np.average(x, axis=0, weights=self.weights)
        self.scale = np.sqrt(np.average((x-self.mean)**2, axis=0, weights=self.weights))
        self.mean[0], self.scale[0] = 0., 1.
        self.scale[self.scale == 0] = 1.
        z = np.array((x-self.mean)/self.scale, order='F')
        root_weight = np.sqrt(self.weights)
        z *= root_weight[:, None]
        q, r = qr(z, mode='economic', overwrite_a=True, check_finite=False)
        u, s, vt = svd(r, full_matrices=False, check_finite=False)
        self.rank = int(np.sum(s > s[0]*rcond))
        self.condition = float(s[0]/s[self.rank-1])
        # Usually full rank: retain QR coordinates to avoid an extra large GEMM.
        if self.rank == x.shape[1]:
            self.q = q
            self.back = lstsq(r, np.eye(len(r)), cond=rcond)[0]
        else:
            self.q = q @ u[:, :self.rank]
            self.back = vt[:self.rank].T/s[:self.rank]
        self.q *= np.sqrt(len(x))
        # q represents the original, unweighted prediction. Base weights belong
        # outside rho; multiplying robust residuals by sqrt(w) changes thresholds.
        self.q /= root_weight[:, None]
        self.back *= np.sqrt(len(x))
        self.n = len(x)

    def native(self, beta):
        scaled = self.back @ beta
        coef = scaled / self.scale
        coef[0] -= self.mean @ coef
        return coef

    def fit(self, target, sigma=None, delta=1.5, max_iter=300):
        target = np.asarray(target, dtype=float)
        if target.shape != (self.n,) or not np.isfinite(target).all():
            raise ValueError('Invalid target')
        beta = self.q.T @ (self.weights*target) / self.n
        if sigma is None:
            return self.native(beta), dict(loss='squared', iterations=0, converged=True,
                                           rank=self.rank, condition=self.condition)
        if not np.isfinite(sigma) or sigma <= 0 or not np.isfinite(delta) or delta <= 0:
            raise ValueError('Positive finite residual scale and Huber threshold required')
        # Centering only improves arithmetic; the constant remains freely fitted.
        center = float(np.average(target, weights=self.weights))
        t = (target-center)/sigma
        start = self.q.T @ (self.weights*t)/self.n
        result = minimize(loss_gradient, start, args=(self.q, t, delta, self.weights), jac=True,
                          method='L-BFGS-B', options=dict(maxiter=max_iter, gtol=1e-8,
                                                        ftol=1e-14, maxls=40))
        grad = float(np.max(np.abs(result.jac)))
        if not result.success or grad > 2e-6:
            raise RuntimeError(f'Huber fit did not converge: {result.message}; gradient={grad:g}')
        coef = self.native(result.x)*sigma
        coef[0] += center
        residual = (self.q @ result.x-t)
        weights = 1/np.hypot(1., residual/delta)
        combined = self.weights*weights
        return coef, dict(loss='pseudo_huber', delta=delta, sigma=float(sigma),
                          iterations=int(result.nit), converged=True, gradient=grad,
                          rank=self.rank, condition=self.condition,
                          fraction_weight_below_half=float(np.mean(weights < .5)),
                          effective_n=float(combined.sum()**2/(combined @ combined)))


def robust_scale(residual):
    residual = np.asarray(residual)
    sigma = 1.482602218505602*np.median(np.abs(residual-np.median(residual)))
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError('Cannot estimate a positive training residual scale')
    return float(sigma)
