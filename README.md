## Gauss-Legendre Quadrature Rule in 1D and 2D

In general a quadrature rule is an approximation of the definite integral of a function `f(x)` over an interval `[a, b]` by a sum of weighted functions at certain points within the domain of integration. The goal is to find such weights and points within the domain. In general, a quadrature is in the form of:

\[\int\limits_{a}^{b} f(x) dx \approx \sum_{i=1}^{n} w_i f(x_i)$\]

where $w_i$'s are the weights and $x_i$'s are the fixed points (the **roots** of a polynomial belonging to a class of orthogonal polynomials) within the domain. In the case of Gaussian quadrature, $a=-1, b=1$, and the quadrature is exact (approximation becomes equality) up to polynomials of degree $2n-1$. Gaussian-Legendre quadrature is a special form of Gaussian quadrature that utilizes orthogonal polynomials that are called Legendre polynomials, denoted by $P_n(x)$. These polynomials are defined as an orthogonal system (orthogonal with respect to $L_2$ inner product) satisfying:

\begin{equation}
<P_m(x), \ P_n(x)> \ = \int\limits_{-1}^{1} P_m(x) P_n(x) dx = 0 \ \ \text{if} \ \ n\neq m
\end{equation}

Some of the first Legendre polynomials are: $P_0(x) = 1, \quad P_1(x) = x, \quad P_2(x) = \frac{3x^2-1}{2} \dots$. Rest of the legendre polynomials can be found using **Bonnet's recursion formula**:

\begin{equation}
(n+1) P_{n+1}(x) = (2n+1)x P_n(x) - n P_{n-1}(x)
\end{equation}


or **Rodrigues' formula**:

\begin{equation}
P_n(x) = \dfrac{1}{2^n n!} \dfrac{d^n}{dx^n} (x^2-1)^n
\end{equation}

With the *n-th* polynomial normalized to given $P_n(1) = 1$, the *i-th* Gauss node $x_i$ is the *i-th* root of $P_n(x)$, and the weights are given by:

\begin{equation}
w_i = \dfrac{2}{(1-x_i^2)({P'}_n(x_i))^2}
\end{equation}

When the integral is taken over a random interval $[a, b]$, interval can easily be changed to $[-1, 1]$ following way:

\begin{equation}
\int\limits_{a}^{b} f(x) dx = \dfrac{b-a}{2} \int\limits_{-1}^{1} f\Big(\frac{b-a}{2}z + \frac{b+a}{2}\Big) dz = c \int\limits_{-1}^{1} f\Big(cz + d\Big) dz
\end{equation}

where $c=\dfrac{b-a}{2}$,  $\ d=\dfrac{b+a}{2}$. Similarly, gaussian quadrature can be generalized to any interval using the transformation above:

\begin{equation}
\int\limits_{a}^{b} f(x) dx = c \int\limits_{-1}^{1} f\Big(cz + d\Big) dz \approx c \sum\limits_{i=1}^{n} w_i f(cx_i +d) = \sum\limits_{i=1}^{n} \tilde{w_i} f(cx_i +d)
\end{equation}

where $\tilde{w_i} = c * w_i$.

### Newton's Method for Root Finding

Newthon-Raphson method is a root-finding algorithm producing successively better approximations to the zeroes (roots) of a real-valued function $f(x)$. Starting with an initial guess $x_0$ (assuming it is close enough), the iteration is given by the formula approximates to a root of the function $f(x)$:

\begin{equation}
x_{n+1} = x_{n} - \dfrac{f(x_n)}{f'(x_n)}
\end{equation}

We will use this method in order to approximate the roots of legendre polynomials and find $x_i$ values for the Gauss-Legendre quadrature.
