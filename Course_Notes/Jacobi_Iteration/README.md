# Jacobi Iteration Method:

For equation set $$\mathbf{Ax} = \mathbf{b}$$, we can decompose $$\mathbf{A} = \mathbf{L} + \mathbf{D} + \mathbf{U}$$

$$\mathbf{L} = \begin{bmatrix}
0 & 0 & 0 & \dots & 0 \\
a_{21} & 0 & 0 & \dots & 0 \\
a_{31} & a_{32} & 0 & \dots & 0 \\
\vdots & \vdots & \vdots &  & \vdots \\
a_{n1} & a_{n2} & a_{n3} & \vdots & 0
\end{bmatrix}$$

$$\mathbf{D} = \begin{bmatrix}
a_{11} &  &  &  &  \\
 & a_{22} &  &  &  \\
 &  & a_{33} &  &  \\
 &  &  & \ddots &  \\
 &  &  &  & a_{nn}
\end{bmatrix}$$

$$\mathbf{U} = \begin{bmatrix}
0 & a_{12} & a_{13} & \dots & a_{1n} \\
0 & 0 & a_{23} & \dots & a_{2n} \\
0 & 0 & 0 & \dots & a_{3n} \\
\vdots & \vdots & \vdots & & \vdots \\
0 & 0 & 0 & \vdots & 0
\end{bmatrix}$$

Then we can get:

$$
\begin{cases}
\mathbf{(L + D + U)x} = \mathbf{b} \\
\mathbf{Dx} = -\mathbf{(L + U)x} + \mathbf{b}\\
\mathbf{x} = -\mathbf{{D}^{-1}(L + U)x} + \mathbf{{D}^{-1}b}
\end{cases}
$$



