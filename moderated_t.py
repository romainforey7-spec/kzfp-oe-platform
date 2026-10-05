"""Empirical-Bayes moderated t-statistics (Smyth, Stat Appl Genet Mol Biol 2004).

A faithful port of the estimators limma uses in `fitFDist`, `squeezeVar` and
`eBayes`. Written because Bioconductor's compiled libraries cannot be loaded in
this sandbox; the accompanying R script (`apms_limma.R`) runs the canonical
implementation and should be used for the archival analysis.

Validated in `validate()` against data simulated from the prior the model
assumes, checking recovery of the prior degrees of freedom and prior variance.
"""
import numpy as np
from scipy.special import digamma, polygamma
from scipy import stats


def trigamma_inverse(x):
    """Solve trigamma(y) = x for y (limma's trigammaInverse: Newton iteration)."""
    x = np.asarray(x, dtype=float)
    y = 0.5 + 1.0 / x                      # limma's starting value
    y = np.where(x > 1e7, 1.0 / np.sqrt(x), y)
    y = np.where(x < 1e-6, 1.0 / x, y)
    for _ in range(50):
        tri = polygamma(1, y)
        dif = tri * (1 - tri / x) / polygamma(2, y)
        y = y + dif
        if np.max(np.abs(dif / y)) < 1e-8:
            break
    return y


def fit_f_dist(var, df1):
    """Estimate the scaled-F prior (d0, s0^2) from per-gene variances."""
    var = np.asarray(var, float)
    df1 = np.broadcast_to(np.asarray(df1, float), var.shape)
    ok = np.isfinite(var) & (var > 0) & (df1 > 0)
    v, d = var[ok], df1[ok]
    z = np.log(v)
    e = z - digamma(d / 2.0) + np.log(d / 2.0)
    ebar = e.mean()
    n = len(e)
    evar = np.sum((e - ebar) ** 2) / (n - 1) - np.mean(polygamma(1, d / 2.0))
    if evar > 0:
        df2 = 2.0 * trigamma_inverse(evar)
        s20 = float(np.exp(ebar + digamma(df2 / 2.0) - np.log(df2 / 2.0)))
    else:
        df2 = np.inf
        s20 = float(np.exp(ebar))
    return float(df2), s20


def squeeze_var(var, df):
    """Shrink per-gene variances toward the fitted prior."""
    df2, s20 = fit_f_dist(var, df)
    if np.isinf(df2):
        post = np.full_like(np.asarray(var, float), s20)
    else:
        post = (df2 * s20 + np.asarray(df, float) * np.asarray(var, float)) / (df2 + np.asarray(df, float))
    return post, df2, s20


def lm_fit(Y, X):
    """Least-squares fit of every row of Y on design X."""
    Y = np.asarray(Y, float); X = np.asarray(X, float)
    n, p = X.shape
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = Y @ X @ XtX_inv.T
    resid = Y - beta @ X.T
    df_res = n - np.linalg.matrix_rank(X)
    sigma2 = (resid ** 2).sum(axis=1) / df_res
    return dict(beta=beta, sigma2=sigma2, df_res=df_res, XtX_inv=XtX_inv, rank=np.linalg.matrix_rank(X))


def contrast_test(fit, contrast):
    """Moderated t-test for one contrast vector over the design columns."""
    c = np.asarray(contrast, float)
    logfc = fit['beta'] @ c
    v = float(c @ fit['XtX_inv'] @ c)                      # unscaled variance factor
    post, df0, s20 = squeeze_var(fit['sigma2'], fit['df_res'])
    se = np.sqrt(post * v)
    t = logfc / se
    df_total = fit['df_res'] + (df0 if np.isfinite(df0) else 1e6)
    p = 2 * stats.t.sf(np.abs(t), df_total)
    return dict(logFC=logfc, t=t, P=p, df_prior=df0, s2_prior=s20,
                df_total=df_total, se=se, sigma2_post=post)


def bh(p):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); ranked = p[o] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[o] = np.clip(ranked, 0, 1)
    return out


def validate(seed=0, n_genes=8000, d_res=7, d0_true=4.0, s20_true=0.25):
    """Simulate from the assumed prior and check the estimators recover it."""
    rng = np.random.default_rng(seed)
    sigma2_g = s20_true * d0_true / rng.chisquare(d0_true, n_genes)
    s2_g = sigma2_g * rng.chisquare(d_res, n_genes) / d_res
    d0_hat, s20_hat = fit_f_dist(s2_g, np.full(n_genes, d_res))
    # null calibration: moderated t should be uniform under the null
    X = np.zeros((d_res + 2, 2)); X[:, 0] = 1; X[(d_res + 2) // 2:, 1] = 1
    Y = rng.normal(0, np.sqrt(sigma2_g)[:, None], size=(n_genes, d_res + 2))
    f = lm_fit(Y, X); r = contrast_test(f, [0, 1])
    ks = stats.kstest(r['P'], 'uniform')
    return dict(d0_true=d0_true, d0_hat=round(d0_hat, 3),
                s20_true=s20_true, s20_hat=round(s20_hat, 4),
                null_p_uniform_KS_stat=round(ks.statistic, 4),
                null_p_uniform_KS_p=round(ks.pvalue, 4),
                frac_p_below_05=round(float((r['P'] < 0.05).mean()), 4))
