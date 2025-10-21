# Assignment 2

## Exercise 3.24

Use Jacobi method, Gauss-Seidel method, and SOR method to solve the following linear system:

$$\begin{bmatrix}
4 & 3 & 0 \\
3 & 4 & -1 \\
0 & -1 & 4
\end{bmatrix}\begin{bmatrix}
x_1 \\ x_2 \\ x_3
\end{bmatrix} = \begin{bmatrix}
24 \\ 30 \\ -24
\end{bmatrix}$$

---

### Tasks:

1. Manual Calculation
   Derive and write out the iteration formulas by hand. Take a photo of your handwritten derivation and include it as an image.

2. Program Implementation and Verification
   - Write code to implement and verify the methods.
   - Include the full source code.
   - Output the result of each iteration step.

1. Test the methods using two different initial guesses and compare the results.
   Restrictions Do not use built-in library functions for solving linear systems. You must implement the core algorithms manually.

---

## Mathematical Deductions

We are given the system:

$$\begin{bmatrix}
4 & 3 & 0 \\
3 & 4 & -1 \\
0 & -1 & 4
\end{bmatrix}\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix} = \begin{bmatrix}
24 \\
30 \\
-24
\end{bmatrix}$$

Denoting:

$$
\mathbf{A} = \begin{bmatrix}
4 & 3 & 0 \\
3 & 4 & -1 \\
0 & -1 & 4
\end{bmatrix},   
\mathbf{x} = \begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix},   
\mathbf{b} = \begin{bmatrix}
24 \\
30 \\
-24
\end{bmatrix}
$$

### Decomposing Matrix $$\( \mathbf{A} \)$$

We know that $$\( \mathbf{A} = \mathbf{L} + \mathbf{D} + \mathbf{U} \)$$, so that:


$$
\mathbf{D} = \begin{bmatrix}
4 & 0 & 0 \\
0 & 4 & 0 \\
0 & 0 & 4
\end{bmatrix}
$$

$$
\mathbf{L} = \begin{bmatrix}
0 & 0 & 0 \\
3 & 0 & 0 \\
0 & -1 & 0
\end{bmatrix}
$$

$$
\mathbf{U} = \begin{bmatrix}
0 & 3 & 0 \\
0 & 0 & -1 \\
0 & 0 & 0
\end{bmatrix}
$$

### 1: Jacobi Iteration Formula

The Jacobi iteration is:

$$
\mathbf{x}^{(k+1)} = \mathbf{D}^{-1} \left( \mathbf{b} - (\mathbf{L} + \mathbf{U}) \mathbf{x}^{(k)} \right)
$$

---

### 2: Gauss-Seidel Iteration Formula

The Gauss-Seidel iteration is:

$$
\mathbf{x}^{(k+1)} = (\mathbf{D} + \mathbf{L})^{-1} \left( \mathbf{b} - \mathbf{U} \mathbf{x}^{(k)} \right)
$$

---

### 3: SOR Iteration Formula

The SOR (Successive Over-Relaxation) iteration is:

$$
\mathbf{x}^{(k+1)} = (1 - \omega)\mathbf{x}^{(k)} + \omega (\mathbf{D} + \mathbf{L})^{-1} \left( \mathbf{b} - \mathbf{U} \mathbf{x}^{(k)} \right)
$$

where $$\( \omega \in (1.8, 1.22) \)$$ is the relaxation factor.
