# MatOpsWrapper Class Documentation

`MatOpsWrapper` is a MATLAB class for performing, tracking, and managing elementary and compound matrix operations. It supports row operations, matrix augmentation, undoing operations, simplification, and LU decomposition, with optional augmentation for solving systems or tracking transformations.

## Constructor

```
obj = MatOpsWrapper(mat)
```
Creates a new wrapper for the matrix `mat`.

## Static builders

- `build(seq, vars)` — Create a `MatOpsWrapper` from a 2-D sequence of symbolic expressions. `vars` is optional; if omitted the method auto-detects symbolic variables. Returns a wrapper whose `mat` contains coefficients for each variable.


## Methods

### Display and State
- `show()`: Displays the current matrix (and augmented matrix if present).
- `steps()`: Displays the list of operations performed on the matrix.
- `freezeMat()`: Returns a copy of the current matrix.
- `freezeOps()`: Returns a copy of the current operations.
- `freezeAug()`: Returns a copy of the current augmentation.
- `freeze()`: Returns a deep copy of the current `MatOpsWrapper` object.

### Elementary Row Operations
- `swapRows(r1, r2)`: Swaps rows `r1` and `r2`.
- `addRow(r1, r2, c)`: Adds `c` times row `r2` to row `r1`.
- `mulRow(r1, c)`: Multiplies row `r1` by scalar `c`.

### Compound Operations
- `ref()`: Reduces the matrix to Row Echelon Form (REF).
- `rref()`: Reduces the matrix to Reduced Row Echelon Form (RREF).

### Operation Management
- `undo(upto)`: Undoes the last `upto` operations.
- `delOps(upto)`: Deletes the first `upto` operations from the operation history.
- `opsMop()`: Returns a `MatOpsWrapper` object of the matrix representing the compound operation.

### Augmentation
- `augment(mat)`: Augments the current matrix with another matrix `mat`.

### Substitution and Simplification
- `subs(match, replacement)`: Substitutes `match` (a symbolic variable) with `replacement` (a numeric value) in the matrix (and augmentation if present).
- `simp()`: Simplifies the matrix (and augmentation if present).

### Matrix Inversion and Decomposition
- `invert()`: Returns the inverse of the matrix as a `MatOpsWrapper` object.
- `lu()`: Returns the LU decomposition as `[l, u]`, where `l` is the `MatOpsWrapper` object of the unit lower triangular matrix and `u` is the `MatOpsWrapper` object of the upper triangular matrix.

### Column Operations
- `subCol(colNum, with)`: Substitutes out a column from the matrix while updating the augmentation.

	- Parameters:
		- `colNum` (integer): Index of the column in the main matrix to remove (1-based).
		- `with` : value to substitute with

	- Returns
		- `MatOpsWrapper` object with the given column substituted out

	- Notes:
		- The method expects the instance to be augmented (i.e., `aug` is a `MatOpsWrapper`) because it updates `aug.mat`. If `aug` is `false`, calling `subCol` will error.

## Example Usage

```matlab
A = [1 2; 3 4];
mop = MatOpsWrapper(A);
mop.swapRows(1, 2);
mop.show();
mop.ref();
mop.show();
```

See `MatOpsWrapper.m` for full implementation details.
