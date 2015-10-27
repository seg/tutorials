import numpy as np
import scipy.ndimage as ndimage
import scipy.signal

def bahorich_coherence(data, zwin):
    ni, nj, nk = data.shape
    out = np.zeros_like(data)
    padded = np.pad(data, ((0, 0), (0, 0), (zwin//2, zwin//2)), mode='reflect')

    for i, j, k in np.ndindex(ni - 1, nj - 1, nk - 1):
        center_trace = data[i,j,:]
        center_std = center_trace.std()
        x_trace = padded[i+1, j, k:k+zwin]
        y_trace = padded[i, j+1, k:k+zwin]

        xcor = np.correlate(center_trace, x_trace)
        ycor = np.correlate(center_trace, y_trace)

        px = xcor.max() / (xcor.size * center_std * x_trace.std())
        py = ycor.max() / (ycor.size * center_std * y_trace.std())
        out[i,j,k] = np.sqrt(px * py)

    return out

def moving_window(data, window, func):
    wrapped = lambda region: func(region.reshape(window))
    return ndimage.generic_filter(data, wrapped, window, mode='reflect')

def marfurt_semblance(region):
    region = region.reshape(-1, region.shape[-1])
    ntraces, nsamples = region.shape

    cov = region.dot(region.T)
    return cov.sum() / cov.diagonal().sum() / ntraces

def semblance2(region):
    region = region.reshape(-1, region.shape[-1])
    ntraces, nsamples = region.shape

    square_of_sums = np.sum(region, axis=0)**2
    sum_of_squares = np.sum(region**2, axis=0)
    return square_of_sums.sum() / sum_of_squares.sum() / ntraces

def eig(region):
    region = region.reshape(-1, region.shape[-1])

    cov = region.dot(region.T)
    vals = np.linalg.eigvalsh(cov)
    return vals.max() / vals.sum()

def complex_semblance(region):
    region = region.reshape(-1, region.shape[-1])
    ntraces, nsamples = region.shape

    region = scipy.signal.hilbert(region)
    region = np.hstack([region.real, region.imag])
    cov = region.dot(region.T)
    return np.abs(cov.sum() / cov.diagonal().sum()) / ntraces

def complex_eig(region):
    region = region.reshape(-1, region.shape[-1])

    region = scipy.signal.hilbert(region)
    region = np.hstack([region.real, region.imag])

    cov = region.dot(region.T)
    vals = np.linalg.eigvals(cov)
    return np.abs(vals.max() / vals.sum())

def flatten(data, surface, window):
    surface = ndimage.gaussian_filter(surface.astype(float), 3)

    ni, nj, nk = data.shape
    ik = np.arange(nk)
    out_ik = np.arange(window) - window // 2

    out = np.zeros((ni, nj, window))
    for i, j in np.ndindex(ni, nj):
        trace = data[i,j,:]
        k = surface[i, j]
        shifted = np.interp(out_ik + k, ik, trace)

        out[i,j,:] = shifted

    return out

def unflatten(data, surface, orig_shape):
    out = np.zeros(orig_shape)
    surface = np.clip(surface, 0, orig_shape[-1] - 1)

    win = data.shape[-1] // 2
    for i, j in np.ndindex(orig_shape[0], orig_shape[1]):
        k = surface[i,j]

        outmin, outmax = max(0, k - win), min(orig_shape[-1], k + win + 1)
        inmin, inmax = outmin - (k - win), k + win + 1 - outmax
        inmax = data.shape[-1] - abs(inmax)

        out[i, j, outmin:outmax] = data[i, j, inmin:inmax]

    return out

def gradients(seismic, sigma):
    """Builds a 4-d array of the gaussian gradient of *seismic*."""
    grads = []
    for axis in range(3):
        # Gaussian filter with order=1 is a gaussian gradient operator
        grad = scipy.ndimage.gaussian_filter1d(seismic, sigma, axis=axis, order=1)
        grads.append(grad[..., np.newaxis])
    return np.concatenate(grads, axis=3)

def moving_window4d(grad, window, func):
    """Applies the given function *func* over a moving *window*, reducing
    the input *grad* array from 4D to 3D."""
    # Pad in the spatial dimensions, but leave the gradient dimension unpadded.
    half_window = [(x // 2, x // 2) for x in window] + [(0, 0)]
    padded = np.pad(grad, half_window, mode='reflect')

    out = np.empty(grad.shape[:3], dtype=float)
    for i, j, k in np.ndindex(out.shape):
        region = padded[i:i+window[0], j:j+window[1], k:k+window[2], :]
        out[i,j,k] = func(region)
    return out

def gst_coherence_calc(region):
    """Calculate gradient structure tensor coherence on a local region.
    Intended to be applied with *moving_window4d*."""
    region = region.reshape(-1, 3)
    gst = region.T.dot(region) # This is the 3x3 gradient structure tensor

    # Reverse sort of eigenvalues of the GST (largest first)
    eigs = np.sort(np.linalg.eigvalsh(gst))[::-1]

    return (eigs[0] - eigs[1]) / (eigs[0] + eigs[1])

def gst_coherence(seismic, window, sigma=1):
    """Randen, et al's (2000) Gradient Structure Tensor based coherence."""
    # 4-d gradient array (ni x nj x nk x 3)
    grad = gradients(seismic, sigma)
    return moving_window4d(grad, window, gst_coherence_calc)
