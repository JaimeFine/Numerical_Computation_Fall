# Assignment 4

Find the smallest positive root of the nonlinear equation:

$$\cos x + \frac{1}{1 + {e}^{-2x}} = 0$$

Take the initial value $$\mathbf{x_{0} = 3}$$, and examine the following iterative schemes:

1. $$x_{k+1} = \arccos \left(-\frac{1}{1 + {e}^{-2x_{k}}} \right)$$ for $$k = 0, 1, 2, ...$$;
2. $$\vcenter x_{k+1} = 0.5\ln \left(-\frac{1}{1+1 \ \cos x_{k}} \right)$$
3. Netwon's method;

For each of the above iterative schemes, first prove that it indeed corresponds to an equivalent fixed-point problem, then theoretically analyze whether it converges locally and what the rate of convergence is. Finally, implement the method and verify your conclusions.

> Remeber to verify the correctness of the solution.