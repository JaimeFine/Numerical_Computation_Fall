# Assignment 3


Randomly generate a $$\mathbf{5000 \times 5000}$$ matrix $$\mathbf{A}$$ and a $$\mathbf{5000 \times 1}$$ vector $$\mathbf{x}$$, then compute $$\mathbf{b} = A\mathbf{x}$$, so that we can solve linear equations of the form $$\mathbf{b =} A\mathbf{x}$$ of any size.

- Solve the equation $$\mathbf{b =} A\mathbf{x}$$ using *library functions*.
- Perform *LU decomposition* on A, then solve for x.
  Compute the norm $$\left\lVert \mathbf{b =} A\mathbf{x} \right\rVert$$ to verify the correctness of the solution.
- Compare the *computation time* between the two methods and analyze the results.
- Then solve the following system:

$$
\mathbf{A} = \begin{bmatrix}
0.8147 & 0.0975 & 0.1576 & 0.1419 & 0.6557 \\
0.9058 & 0.2785 & 0.9706 & 0.4218 & 0.0357 \\
0.1270 \cdot 10^{10} & 0.5469 & 0.9572 & 0.9157 & 0.8491 \\
0.9134 & 0.9575 & 0.4854 \cdot 10^{8} & 0.7922 & 0.9340 \\
0.6324 & 0.9649 & 0.8003 & 0.9595 & 0.6787
\end{bmatrix}
$$

$$
\mathbf{b} = \mathbf{1.0}^{9} \cdot \begin{bmatrix}
0.000000002258000 \\
0.000000001597700 \\
1.270000002354900 \\
0.024270003904200 \\
0.000000003360250
\end{bmatrix}
$$

Find the corresponding solution x.

> **Remember to verify the correctness of the solution by computing the residual norm!**