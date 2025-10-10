# Newton-Raphson Method:

**$$x_{\text{new}} = x - \frac{f(x)}{f'(x)}$$**

## The algorithm:
- Functionality: Find the root of the function f(x) = 0
- Input: f(x), f'(x), x0, max_iter, e
- Output: The approximate root, or failed message

1. flag <- 0, count <- 0, e <- (1 - |f'(x0)|) * e
2. x1 <- x0 - f(x0) / f'(x0); count <- count + 1
3. if |x1 - x0| <= e, then flag <- 1, JMP 4;
   else if count <= max_iter, then x0 <- x1, JMP 2;
   else, JMP 4
4. if flag = 1, then c <- x1, print c;
   else, print failed message

## A simplified algorithm:
- Functionality: Find the root of the function f(x) = 0
- Input: f(x), f'(x), x0, max_iter, e
- Output: The approximate root, or failed message

1. flag <- 0, count <- 0, e <- (1 - |f'(x0)|) * e, x <- x0
2. x1 <- x0 - f(x0) / f'(x); count <- count + 1
3. if |x1 - x0| <= e, then flag <- 1, JMP 4;
   else if count <= max_iter, then x0 <- x1, JMP 2;
   else, JMP 4
4. if flag = 1, then c <- x1, print c;
   else, print failed message