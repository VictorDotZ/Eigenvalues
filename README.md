# Eigenvalues
QR decomposition using Givens rotations and finding eigenvalues using Householder Reflections

## Compilation
```bash
g++ main.cpp -O0
```

## Usage
```bash
./a.out <matrix size> <output size> <precision> <formula>
```

for example:

```bash
./a.out 500 5 1e-10 3
```

---

> Tip: compile with -O0 to prevent compiler optimizations and to demonstrate the asymptotic behavior of the algorithm with smaller matrix sizes.

> Tip: use less precision for first and second formulas.
