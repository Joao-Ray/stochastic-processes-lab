"""Independent likelihood, filtering, Fisher-information and KL references."""
import numpy as np
import pytest
from scipy.linalg import null_space
from scipy.stats import multivariate_normal
from src.ou_measurement_error import (prepare_sample, profile_shape, fit_noisy_ou,
    profile_phi, profile_membership, information_matrix, information_limit)


@pytest.mark.parametrize('repeats,gamma', [(1, 0.), (1, .6), (1, 1.), (2, .6), (4, .9)])
@pytest.mark.parametrize('restricted', [False, True])
def test_profile_against_dense_gaussian_and_contrasts(repeats, gamma, restricted):
    rng = np.random.default_rng(42)
    sites, phi = 8, .87
    y = 2+rng.normal(size=(sites, repeats))
    distance = np.abs(np.arange(sites)[:, None]-np.arange(sites))
    covariance = gamma*np.kron(phi**distance, np.ones((repeats, repeats)))+(1-gamma)*np.eye(y.size)
    inverse = np.linalg.inv(covariance)
    ones = np.ones(y.size)
    mu = ones @ inverse @ y.ravel()/(ones @ inverse @ ones)
    residual = y.ravel()-mu
    degrees = y.size-int(restricted)
    scale = residual @ inverse @ residual/degrees
    expected = .5*(np.linalg.slogdet(covariance)[1]+degrees*(np.log(2*np.pi*scale)+1))
    if restricted:
        expected += .5*np.log(ones @ inverse @ ones)
        contrasts = null_space(ones[None, :])
        c = contrasts.T @ covariance @ contrasts
        projected = contrasts.T @ y.ravel()
        contrast_nll = .5*(np.linalg.slogdet(c)[1]+degrees*(np.log(2*np.pi*(projected @ np.linalg.solve(c, projected))/degrees)+1))
        assert expected == pytest.approx(contrast_nll+.5*np.log(y.size), abs=1e-10)
    actual = profile_shape(y, phi, gamma, restricted)
    assert actual['nll'] == pytest.approx(expected, abs=1e-9)
    assert actual['mu'] == pytest.approx(mu, abs=1e-10)
    assert actual['scale'] == pytest.approx(scale, abs=1e-10)


@pytest.mark.parametrize('repeats', [1, 3])
def test_profile_ml_against_sequential_kalman_likelihood(repeats):
    y = np.random.default_rng(7).normal(size=(13, repeats))
    phi, gamma = .94, .65
    fit = profile_shape(y, phi, gamma)
    latent, noise = gamma*fit['scale'], (1-gamma)*fit['scale']
    mean, var, nll = fit['mu'], latent, 0.
    for site in y:
        for value in site:
            error, innovation = value-mean, var+noise
            nll += .5*(np.log(2*np.pi*innovation)+error**2/innovation)
            gain = var/innovation
            mean += gain*error
            var *= 1-gain
        mean = fit['mu']+phi*(mean-fit['mu'])
        var = phi*phi*var+latent*(1-phi*phi)
    assert fit['nll'] == pytest.approx(nll, abs=1e-9)


def test_fisher_against_full_replicated_covariance_derivatives():
    theta, v, eta2, dt, m, k = .8, 1.2, .7, .2, 6, 3
    distance = np.abs(np.arange(m)[:, None]-np.arange(m))*dt
    def covariance(par):
        t, variance, noise = par
        return variance*np.kron(np.exp(-t*distance), np.ones((k, k)))+noise*np.eye(m*k)
    par = np.array([theta, v, eta2])
    base = covariance(par)
    derivatives = []
    for j in range(3):
        h = np.zeros(3); h[j] = 1e-5
        derivatives.append((covariance(par+h)-covariance(par-h))/(2e-5))
    products = [np.linalg.solve(base, d) for d in derivatives]
    expected = np.array([[.5*np.trace(a @ b) for b in products] for a in products])
    actual = information_matrix(theta, v, eta2, dt, m, k)
    np.testing.assert_allclose(actual['matrix'], expected, atol=2e-8, rtol=2e-8)
    assert actual['efficient_theta'] > 0


def test_kl_bound_against_direct_gaussian_formula_and_white_limit():
    m, phi0, phi1, gamma = 12, .9, .7, .2
    d = np.abs(np.arange(m)[:, None]-np.arange(m))
    p, q = gamma*phi0**d+(1-gamma)*np.eye(m), gamma*phi1**d+(1-gamma)*np.eye(m)
    expected = .5*(np.trace(np.linalg.solve(q, p))-m+np.linalg.slogdet(q)[1]-np.linalg.slogdet(p)[1])
    result = information_limit(phi0, phi1, gamma, m)
    assert result['kl_0_to_1'] == pytest.approx(expected, abs=1e-12)
    assert result['analytic_kl_upper'] >= expected
    white = information_limit(phi0, phi1, 0., m)
    assert white['tv_upper'] < 1e-12
    assert white['simultaneous_membership_lower'] == pytest.approx(.9)
    assert white['worst_error_probability_lower'] == pytest.approx(.5)


def test_affine_equivariance_of_shape_and_fit():
    y = np.random.default_rng(8).normal(size=(20, 2))
    for restricted in [False, True]:
        first = profile_shape(y, .8, .6, restricted)
        shifted = profile_shape(7-3*y, .8, .6, restricted)
        assert shifted['nll']-first['nll'] == pytest.approx((y.size-int(restricted))*np.log(3), abs=1e-9)
        assert shifted['mu'] == pytest.approx(7-3*first['mu'], abs=1e-10)
        assert shifted['scale'] == pytest.approx(9*first['scale'], abs=1e-10)
    fit = fit_noisy_ou(y, .2)
    assert fit['white_limit'] or fit['search_cap'] or fit['identified'] or not fit['optimizer_converged']
    assert profile_phi(y, fit['phi'])['nll'] <= fit['nll']+1e-5
    result = profile_membership(y, .2, 1., fit)
    assert result['available'] and result['lr'] >= 0


def test_invalid_samples_and_singular_replicates():
    for y in [[1.]*5, [1., np.nan, 2., 3.], np.ones((5, 2))]:
        with pytest.raises(ValueError): prepare_sample(y)
    with pytest.raises(ValueError): profile_shape([1.,2.,.2,.3], 1., .5)
    with pytest.raises(ValueError): fit_noisy_ou([1.,2.,.2,.3], 0.)
    with pytest.raises(ValueError): information_matrix(1.,1.,0.,.1,8)


def test_white_limit_has_no_invented_drift_or_variance_decomposition():
    fit=fit_noisy_ou([.1,-.6,.3,.8,-.4,.1],.2,fixed_gamma=0.)
    assert fit['white_limit'] and not fit['identified']
    assert fit['theta'] is None and fit['noise_variance'] is None and fit['latent_variance'] is None


def test_noisy_study_cli_replay_and_all_denominators(tmp_path):
    import csv,gzip,hashlib,json,os,subprocess,sys
    from pathlib import Path
    root=Path(__file__).resolve().parents[1]
    folder=root/'subprojects/ou-inference-study/measurement-error'
    config=json.loads((folder/'ci-config.json').read_text())
    config.update(replications=4,audit_records_per_scenario=0)
    path=tmp_path/'config.json';path.write_text(json.dumps(config))
    outputs=[tmp_path/'one',tmp_path/'two']
    for dest in outputs:
        subprocess.run([sys.executable,str(folder/'run.py'),'--config',str(path),'--output',str(dest)],
            cwd=tmp_path,env={**os.environ,'MPLCONFIGDIR':str(tmp_path/'mpl')},check=True,capture_output=True,text=True,timeout=120)
    archive=outputs[0]/'replications.csv.gz'
    assert archive.read_bytes()==(outputs[1]/archive.name).read_bytes()
    summary=json.loads((outputs[0]/'summary.json').read_text())
    assert summary['record_archive_sha256']==hashlib.sha256(archive.read_bytes()).hexdigest()
    for name,digest in summary['source_sha256'].items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest
    with gzip.open(archive,'rt') as f:rows=list(csv.DictReader(f))
    assert len(rows)==16
    for s in summary['scenarios']:
        selected=[r for r in rows if r['scenario']==s['id']]
        assert len(selected)==s['attempted']==4
        for method,metric in s['methods'].items():
            available=[r for r in selected if r[method+'_test_available']=='True']
            hits=sum(r[method+'_hit']=='True' for r in selected)
            assert len(available)==metric['test_available']
            assert metric['success_and_coverage']['value']==hits/4
            if available: assert metric['conditional_coverage']['value']==hits/len(available)
            assert sum(r[method+'_theta']!='' for r in selected)==metric['point']['available']
    assert (outputs[0]/'coverage_and_error.pdf').read_bytes().startswith(b'%PDF')


def test_reml_ridge_regression_matches_independent_search_cap_reference():
    # Record exposed by differential-evolution audit, not a chosen narrow CI.
    from src.ou_monte_carlo import stationary_null_paths
    dt=10/99
    latent=stationary_null_paths(float(np.exp(-dt)),99,500,
        np.random.default_rng(np.random.SeedSequence([20261015,300,4,0])))
    sensor=np.random.default_rng(np.random.SeedSequence([20261015,300,4,1])).standard_normal((500,100,2))
    y=latent[0,:,None]+sensor[0]
    fit=fit_noisy_ou(y,dt,restricted=True)
    assert fit['nll'] <= 291.43653156904224+1e-4
    assert fit['search_cap'] and fit['theta'] is None


def test_nuisance_profiled_starts_recover_missed_low_signal_mode():
    import json
    from pathlib import Path
    root=Path(__file__).resolve().parents[1]
    fixture=json.loads((root/'subprojects/ou-inference-study/measurement-error/development/initial_extreme_profiles.json').read_text())['examples'][0]
    y=np.array(fixture['raw_observations'])
    fit=fit_noisy_ou(y,fixture['dt'],restricted=True)
    moderate=profile_phi(y,float(np.exp(-4*fixture['dt'])),restricted=True)
    assert fit['nll'] < 323.24651292345254-.1
    assert fit['nll'] <= moderate['nll']+1e-6
    assert fit['theta'] is not None and 3 < fit['theta'] < 6
