# %%
from collections.abc import Sequence

def simpson_uniform(x: Sequence[float], f: Sequence[float]) -> float:
    """
    Simpson rule for regularly spaced data.

    :param x: Sampling points for the function values
    :param f: Function values at the sampling points

    :return: approximation for the integral
    """
    N = len(x) - 1
    h = x[1] - x[0]

    result = (h/3) * (f[0] + 2*sum(f[:N-1:2]) + 4*sum(f[1:N:2]) + f[N])
    return result

def simpson_nonuniform(x: Sequence[float], f: Sequence[float]) -> float:
    """
    Simpson rule for irregularly spaced data.

    :param x: Sampling points for the function values
    :param f: Function values at the sampling points

    :return: approximation for the integral
    """
    N = len(x) - 1
    h = [x[i + 1] - x[i] for i in range(0, N)]
    assert N > 0

    result = 0.0
    for i in range(1, N, 2):
        h0, h1 = h[i - 1], h[i]
        hph, hdh, hmh = h1 + h0, h1 / h0, h1 * h0
        result += (hph / 6) * (
            (2 - hdh) * f[i - 1] + (hph**2 / hmh) * f[i] + (2 - 1 / hdh) * f[i + 1]
        )

    print(f"Pre-correction: {result}")

    correction = 0.0
    if N % 2 == 1:
        h0, h1 = h[N - 2], h[N - 1]

        correction += f[N]     * (2 * h1 ** 2 + 3 * h0 * h1) / (6 * (h0 + h1))
        correction += f[N - 1] * (h1 ** 2 + 3 * h1 * h0)     / (6 * h0)
        correction -= f[N - 2] * h1 ** 3                     / (6 * h0 * (h0 + h1))
        result += correction

    print(f"Correction: {correction}")
    print(f"Post-correction: {result}")

    return result

# %%
import numpy as np

a = 0
b = np.pi
n = 12
x = np.linspace(a, b, n)
# x = np.array([0.0, 0.1, 0.5, 0.8, 0.75, 0.9])
f = np.sin(x)

I_uniform = simpson_uniform(x, f)
I_nonuniform = simpson_nonuniform(x, f)

print(I_uniform, I_nonuniform)
# %%
from scipy.integrate import simpson
from scipy.integrate._quadrature import tupleset, _basic_simpson
import numpy as np
# %%
n = 16
a = 0
b = np.pi/2
h = (b - a) / (n - 1)
x = np.linspace(a, b, n)

y = np.sin(x)

simpson(y, x, even='first')
# %%
# %%
axis = -1
dx = 1.0
y = np.asarray(y)
nd = len(y.shape)
N = y.shape[axis]
last_dx = dx
first_dx = dx
returnshape = 0

x = np.asarray(x)
if len(x.shape) == 1:
    shapex = [1] * nd
    shapex[axis] = x.shape[0]
    saveshape = x.shape
    returnshape = 1
    x = x.reshape(tuple(shapex))

val = 0.0
result = 0.0
slice1 = (slice(None),)*nd
slice2 = (slice(None),)*nd

slice1 = tupleset(slice1, axis, -1)
slice2 = tupleset(slice2, axis, -2)
last_dx = x[slice1] - x[slice2]

slice1 = tupleset(slice1, axis, -1)
slice2 = tupleset(slice2, axis, -2)
if x is not None:
    last_dx = x[slice1] - x[slice2]
val += 0.5*last_dx*(y[slice1]+y[slice2])
result = _basic_simpson(y, 0, N-3, x, dx, axis)
# %%
