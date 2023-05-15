abstract type QuadratureMethod end

function integration_weights(mesh::RadialMesh, method::QuadratureMethod)
    weights = similar(mesh)
    return integration_weights!(weights, mesh, method)
end


@doc raw"""
Trapezoidal rule quadrature implemented as described in the
[Wikipedia article](https://en.wikipedia.org/wiki/Trapezoidal_rule).

Non-uniform grid
```math
\int_a^b f(x) dx \approx
\sum_{i=0}^{N} \frac{f(x_{i-1}) + f(x_i)}{2} (x_i - x_{i-1})
```

Uniform grid
```math
\int_a^b f(x) dx \approx
\Delta x \left( \sum_{i=1}^{N-1} f(x_i) + \frac{f(x_N) + f(x_0)}{2} \right)
```
"""
struct Trapezoid <: QuadratureMethod end

function integration_weights!(weights::AbstractVector, mesh::UniformMesh, ::Trapezoid)
    Δx = mesh.a

    weights[begin] = Δx / 2
    weights[(begin + 1):(end - 1)] .= Δx
    weights[end] = Δx / 2

    return weights
end

function integration_weights!(weights::AbstractVector, mesh::RadialMesh, ::Trapezoid)
    weights[begin] = diff(mesh, 1) / 2
    for i in (firstindex(weights) + 1):(lastindex(weights) - 1)
        # Δx[i] + Δx[i-1] = (x[i+1] - x[i]) + (x[i] - x[i-1])
        #                 = x[i+i] - x[i-1]
        weights[i] = (mesh[i + 1] - mesh[i - 1]) / 2
    end
    weights[end] = diff(mesh, mesh.n - 1) / 2

    return weights
end


@doc raw"""
Composite Simpson's rule quadrature implemented as described in the
[Wikipedia article](https://en.wikipedia.org/wiki/Simpson%27s_rule).

## Non-uniform grid

```math
\int_a^b f(x) dx \approx
\sum_{i=0}{N/2-1} \frac{\Delta x_{2i} + \Delta x_{2i+1}}{6} \left[
    \left( 2 - \frac{\Delta x_{2i+1}}{\Delta x_{2i}} \right) f_{2i} +
    \frac{
        \left( \Delta x_{2i} + \Delta x_{2i+1} \right)^2
        }{
        \Delta x_{2i} \Delta x_{2i+1}
    } f_{2i+1} +
    \left( 2 - \frac{\Delta x_{2i}}{\Delta x_{2i + 1}} \right) f_{2i + 2}
\right]
```
where
```math
f_{k} = f \left( a + \sum_{i=0}^{k-1} \Delta x_{i} \right)
```

In the case of an odd number of subintervals, the above formulae are used up to the second
to last interval, and the last interval is computed seperately as
```math
\frac{
    2 \Delta x_{N-1}^2 + 3 \Delta x_{N-1} \Delta x_{N-2}
    }{
    6\left( \Delta x_{N-1} + \Delta x{N-1} \right)
} f_{N} +
\frac{
    \Delta x_{N-1}^2 + 3 \Delta x_{N-1} \Delta x_{N-1}
    }{
    6 \Delta x_{N-1}
} f_{N-1} +
\frac{
    \Delta x_{N-1}^3
    }{
    6 \Delta x_{N-2} \left( \Delta x_{N-2} + \Delta x_{N-1} \right)
} f_{N-2}
```

## Uniform grid

```math
\int_a^b f(x) dx \approx
\frac{1}{3} \Delta x \left[
    f(x_0) + 4 \sum_{i=0}{N/2} f(x_{2i-1}) + 2 \sum_{i=0}{N/2-1} f(x_{2i}) + f(x_N)
\right]
```

In the case of an odd number of subintervals, the above formula is used for the first N-1
subintervals, and the last subinterval is handled with the Trapezoid method. This approach
is generally more well behaved in the specific case of pseudopotential quantities because
the function is generally close to zero in the last interval and the error made by
the Trapezoid rule w.r.t. Simpson's rule has a smaller effect on the value of the integral.
This approach is also equivalent to the approach taken in the non-uniform case.
"""
struct Simpson <: QuadratureMethod end

function integration_weights!(weights::AbstractVector, mesh::UniformMesh, ::Simpson)
    Δx = mesh.a
    N = length(weights) - 1  # Number of intervals

    if !isodd(N)  # Standard Simpson's composite 1/3 rule
        weights[begin] = 1 / 3 * Δx
        weights[(begin + 1):2:(end - 1)] .= 4 * 1 / 3 * Δx
        weights[(begin + 2):2:(end - 1)] .= 2 * 1 / 3 * Δx
        weights[end] = 1 / 3 * Δx
    else # * For functions which decay to zero as r -> ∞, this is a better approximation
        # If the number of intervals is odd, apply Simpsons method to the first N-1 intervals
        # and the Trapezoidal method to the last interval.
        weights[begin] = 1 / 3 * Δx
        weights[(begin+1):2:(end-2)] .= 4 * 1 / 3 * Δx
        weights[(begin+2):2:(end-2)] .= 2 * 1 / 3 * Δx
        weights[end-1] = 5 / 6 * Δx
        weights[end] = 1 / 2 * Δx
    end
    # else
    #     # If the number of intervals is odd, apply Simpsons method to the last N-1 intervals
    #     # and the Trapezoidal method to the first interval.
    #     weights[begin] = 1 / 2 * Δx
    #     weights[begin+1] = 5 / 6 * Δx
    #     weights[(begin+2):2:(end-1)] .= 4 * 1 / 3 * Δx
    #     weights[(begin+3):2:(end-1)] .= 2 * 1 / 3 * Δx
    #     weights[end] = 1 / 3 * Δx
    # end
    # else
    #     # If the number of intervals is odd, average the results of applying Simpson's
    #     # composite 1/3 rule to the first N-1 and last N-1 intervals, applying the
    #     # Trapezoidal rule to the last / first interval.
    #     weights[begin] = 5 / 12 * Δx  # (1/3 + 1/2) / 2 = 5/12
    #     weights[begin + 1] = 13 / 12 * Δx  # (4/3 + 1/3 + 1/2) / 2 = 13/12
    #     weights[(begin + 2):(end - 2)] .= 1 * Δx  # (4/3 + 2/3) / 2 = 1
    #     weights[end - 1] = 13 / 12 * Δx  # ''
    #     weights[end] = 5 / 12 * Δx  # ''
    # end

    return weights
end

function integration_weights!(weights::AbstractVector, mesh::RadialMesh, ::Simpson)
    N = length(weights) - 1  # Number of intervals
    fill!(weights, 0)

    # Skip the last interval if the number of intervals is odd
    istop = isodd(N) ? lastindex(weights) - 3 : lastindex(weights) - 2

    for i in firstindex(weights):2:istop
        Δx_0 = diff(mesh, i)
        Δx_1 = diff(mesh, i + 1)
        prefac = (Δx_0 + Δx_1) / 6
        weights[i] += prefac * (2 - Δx_1 / Δx_0)
        weights[i + 1] += prefac * (Δx_0 + Δx_1)^2 / (Δx_0 * Δx_1)
        weights[i + 2] += prefac * (2 - Δx_0 / Δx_1)
    end

    if isodd(N)  # This handles the last interval when the number of intervals is odd
        Δx_n = diff(mesh, mesh.n - 1)
        Δx_nm1 = diff(mesh, mesh.n - 2)
        weights[end] += (2 * Δx_n^2 + 3 * Δx_n * Δx_nm1) / (6 * (Δx_nm1 + Δx_n))
        weights[end - 1] += (Δx_n^2 + 3 * Δx_n * Δx_nm1) / (6 * Δx_nm1)
        weights[end - 2] -= Δx_n^3 / (6 * Δx_nm1 * (Δx_nm1 + Δx_n))
    end

    return weights
end


@doc raw"""
QuantumESPRESSO Simpson's (1/3) rule quadrature. The expression is equivalent to the
uniform grid case of `Simpson` for an even number of subintervals. However, the
derivative of the mesh `dx` is used in place of the finite difference between adjacent
mesh points `Δx`.
"""
struct QESimpson <: QuadratureMethod end

function integration_weights!(weights::AbstractVector, mesh::RadialMesh, ::QESimpson)
    i = (firstindex(weights) + 1):(lastindex(weights) - 1)
    weights[i] .= 2 / 3 * abs.(mod.(i, 2) .- 2) .* deriv.(Ref(mesh), i)
    if mod(length(weights), 2) == 1
        weights[begin] = 1 / 3 * deriv(mesh, 1)
        weights[end] = 1 / 3 * deriv(mesh, mesh.n - 1)
    else
        weights[begin] = 1 / 3 * deriv(mesh, 1)
        weights[end - 1] = 1 / 3 * deriv(mesh, mesh.n - 2)
    end
    return weights
end


@doc raw"""
ABINIT quadrature -- a collection of composite closed Newton-Cotes rules. Behavior depends
on the number of available points. This method _only_ supports quantities on uniform
meshes!

```math
\int_a^b f(x) dx \approx \begin{cases}
    \frac{\Delta x}{72} \left[
        23.75 f_1 + 95.10 f_2 + 55.20 f_3 + 79.30 f_4 + 70.65 f_5 +
        72 \sum_{i=6}^{N-5} f_i +
        70.65 f_{N-4} 79.30 f_{N-3} + 55.20 f_{N-2} + 95.10 f_{N-1} + 23.75 f_{N}
    \right] & N >= 10 \\
    \frac{\Delta x}{48} \left[
        17 f_1 + 59 f_2 + 43 f_3 + 49 f_4 +
        48 f_5 +
        49 f_6 + 43 f_7 + 59 f_8 + 17 f_9
    \right] & N = 9 \\
    \frac{\Delta x}{48} \left[
        17 f_1 + 59 f_2 + 43 f_3 + 49 f_4 +
        49 f_5 + 43 f_6 + 59 f_7 + 17 f_8
    \right] & N = 8 \\
    \frac{\Delta x}{48} \left[
        17 f_1 + 59 f_2 + 43 f_3 + 50 f_4 + 43 f_5 + 59 f_6 + 17 f_7
    \right] & N = 7 \\
    \frac{\Delta x}{48} \left[
        17 f_1 + 59 f_2 + 44 f_3 + 44 f_4 + 59 f_5 + 17 f_6
    \right] & N = 6 \\
    \frac{\Delta x}{3} \left[
        f_1 + 4 f_2 + 2 f_3 + 4 f_4 + f_5
    \right] & N = 5 \\
    \frac{\Delta x}{8} \left[
        3 f_1 + 9 f_2 + 9 f_3 + 3 f_4
    \right] & N = 4 \\
    \frac{\Delta x}{3} \left[
        1 f_1 + 8 f_2 + 1 f_3
    \right] & N = 3 \\
    \frac{\Delta x}{2} \left[
        f_1 + f_2
    \right] & N = 2 \\
    \frac{\Delta x} \left[
        f_1
    \right] & N = 1
\end{cases}
```
"""
struct AbinitCorrectedTrapezoid <: QuadratureMethod end

function integration_weights!(weights::AbstractVector, mesh::UniformMesh, ::AbinitCorrectedTrapezoid)
    Δx = mesh.a

    if length(weights) >= 10
        weights[begin] = 23.75 / 72 * Δx
        weights[begin + 1] = 95.10 / 72 * Δx
        weights[begin + 2] = 55.20 / 72 * Δx
        weights[begin + 3] = 79.30 / 72 * Δx
        weights[begin + 4] = 70.65 / 72 * Δx
        weights[(begin + 5):(end - 5)] .= Δx
        weights[end - 4] = 70.65 / 72 * Δx
        weights[end - 3] = 79.30 / 72 * Δx
        weights[end - 2] = 55.20 / 72 * Δx
        weights[end - 1] = 95.10 / 72 * Δx
        weights[end] = 23.75 / 72 * Δx
    elseif length(weights) == 9
        weights[begin] = 17 / 48 * Δx
        weights[begin + 1] = 59 / 48 * Δx
        weights[begin + 2] = 43 / 48 * Δx
        weights[begin + 3] = 49 / 48 * Δx
        weights[begin + 4] = Δx
        weights[end - 3] = 49 / 48 * Δx
        weights[end - 2] = 43 / 48 * Δx
        weights[end - 1] = 59 / 48 * Δx
        weights[end] = 17 / 48 * Δx
    elseif length(weights) == 8
        weights[begin] = 17 / 48 * Δx
        weights[begin + 1] = 59 / 48 * Δx
        weights[begin + 2] = 43 / 48 * Δx
        weights[begin + 3] = 49 / 48 * Δx
        weights[end - 3] = 49 / 48 * Δx
        weights[end - 2] = 43 / 48 * Δx
        weights[end - 1] = 59 / 48 * Δx
        weights[end] = 17 / 48 * Δx
    elseif length(weights) == 7
        weights[begin] = 17 / 48 * Δx
        weights[begin + 1] = 59 / 48 * Δx
        weights[begin + 2] = 43 / 48 * Δx
        weights[begin + 3] = 50 / 48 * Δx
        weights[end - 2] = 43 / 48 * Δx
        weights[end - 1] = 59 / 48 * Δx
        weights[end] = 17 / 48 * Δx
    elseif length(weights) == 6
        weights[begin] = 17 / 48 * Δx
        weights[begin + 1] = 59 / 48 * Δx
        weights[begin + 2] = 44 / 48 * Δx
        weights[end - 2] = 44 / 48 * Δx
        weights[end - 1] = 59 / 48 * Δx
        weights[end] = 17 / 48 * Δx
    elseif length(weights) == 5
        weights[begin] = 1 / 3 * Δx
        weights[begin + 1] = 4 / 3 * Δx
        weights[begin + 2] = 2 / 3 * Δx
        weights[end - 1] = 4 / 3 * Δx
        weights[end] = 1 / 3 * Δx
    elseif length(weights) == 4
        weights[begin] = 3 / 8 * Δx
        weights[begin + 1] = 9 / 8 * Δx
        weights[end - 1] = 9 / 8 * Δx
        weights[end] = 3 / 8 * Δx
    elseif length(weights) == 3
        weights[begin] = 1 / 3 * Δx
        weights[begin + 1] = 8 / 3 * Δx
        weights[end] = 1 / 3 * Δx
    elseif length(weights) == 2
        weights[begin] = 1 / 2 * Δx
        weights[end] = 1 / 2 * Δx
    elseif length(weights) == 1
        weights[begin] = Δx
    end
    return weights
end
