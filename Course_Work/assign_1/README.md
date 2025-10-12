# Assignment 1

### Exercise 2.3 — Newton's Method Applications

Use Newton's method to find:

1. The **positive root** of the equation:  
   $$
   x^3 - x^2 - x - 1 = 0
   $$

2. The **smallest positive root** of the equation:  
   $$
   \cos x = \frac{1}{2} + \sin x
   $$

### Convergence Requirement

The iterative difference must satisfy:  
$$
|x_k - x_{k-1}| < \frac{1}{2} \times 10^{-5}
$$

---

Write a program to implement Exercise 2.3. Determine the convergence of the solution to Exercise 2.3, handwrite and photograph the steps, and paste them into this report (only the derivation of the iteration formula is needed). Write a program to verify it; do not call existing library functions, but complete the core algorithm yourself. You do not need to draw a flowchart or write a principle overview. The derivation part should be slightly detailed, and comments should be added to the code to explain parameter settings. The data results should be detailed and aesthetically pleasing; necessary identification and explanatory analysis should be provided for the experimental results.

---

## Mathematic Deductions

### Part 1: Equation  
$$
x^3 - x^2 - x - 1 = 0
$$

Let  
$$
f(x) = x^3 - x^2 - x - 1
$$  
Then  
$$
f(x) \in C^2[0, 2]
$$

#### Root Existence

- $$
  f(1) = 1 - 1 - 1 - 1 = -2 < 0
  $$
- $$
  f(2) = 8 - 4 - 2 - 1 = 1 > 0
  $$

So by the Intermediate Value Theorem, there exists a positive root  
$$
x^* \in [1, 2]
$$  
(The student used the interval \([0, 2]\).)

Also,  
$$
f'(x^*) \ne 0
$$

---

#### Derivatives

- First derivative:  
  $$
  f'(x) = 3x^2 - 2x - 1
  $$

- Second derivative:  
  $$
  f''(x) = 6x - 2
  $$

On the interval \([0, 2]\):  
$$
f''(x) = 6x - 2 \ge 0 \quad \text{for} \quad x \ge \frac{1}{3}
$$

---

#### Newton's Iterative Formula

$$
x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)} = x_k - \frac{x_k^3 - x_k^2 - x_k - 1}{3x_k^2 - 2x_k - 1}
$$

---

### Part 2: Equation  
$$
\cos x = \frac{1}{2} + \sin x
$$

Rewriting:  
$$
\sin x - \cos x + \frac{1}{2} = 0
$$

Let  
$$
f(x) = \sin x - \cos x + \frac{1}{2}
$$

---

#### Trigonometric Identity

Using:  
$$
\sin x - \cos x = \sqrt{2} \sin\left(x - \frac{\pi}{4}\right)
$$

So:  
$$
f(x) = \sqrt{2} \sin\left(x - \frac{\pi}{4}\right) + \frac{1}{2}
$$

---

#### Derivative

$$
f'(x) = \sqrt{2} \cos\left(x - \frac{\pi}{4}\right)
$$

---

#### Root Interval

We seek the **smallest positive root**.  
The student's notes suggest searching in the interval  
$$
[0, \frac{\pi}{4}]
$$

On this interval:  
$$
f'(x) = \sqrt{2} \cos\left(x - \frac{\pi}{4}\right) \ne 0
$$

---

#### Newton's Iterative Formula

$$
x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)} = x_k - \frac{\sqrt{2} \sin\left(x_k - \frac{\pi}{4}\right) + \frac{1}{2}}{\sqrt{2} \cos\left(x_k - \frac{\pi}{4}\right)}
$$
