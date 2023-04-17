# Numerial solution of the Laplace Equation by Jacobi method
### Background
Please refer to [background](./aux/background.md).

---

### Parallelization
Please refer to [parallelization](./aux/parallel.md).

---

**Dependencies**
- Parallel version requires [**mpi compliant library**](https://www.open-mpi.org/)
- Accelerated version requires [**openacc**](https://www.openacc.org/)
- Plotting requires [**gnuplot**](http://www.gnuplot.info/)
- Creating gif requires [**imagemagick**](https://imagemagick.org/image/wizard.png)

## Compilation
Use:

```make [mode]```

with `[mode]`
- Left blank to compile serial version, will produce *jacobi.x* executable
- `mpi` to compile parallel version, will produce *mpi_jacobi.x* executable, requires `mpicc`
- `openacc` to compile parallel accelerated version, will produce *openacc_jacobi.x* executable, requires `pgcc`

## Execution
Use:
```
make run [dim=%d] [iters=%d]
```
to run the serial version, requires `jacobi.x` to be present. 

Use:
```
make mpirun [dim=%d] [iters=%d] [prc=%d]
```
to run the serial version, requires `mpi_jacobi.x` to be present. 

Use:
```
make openaccrun [dim=%d] [iters=%d] [prc=%d]
```
to run the serial version, requires `openacc_jacobi.x` to be present. 

Where:
- `dim` dimension of the internal grid
- `iters` number of iterations
- `prc` number of processes

All will produce *data/solution.dat*.

## Plot
Use:
- `make plot` to see the solution

<img src="./aux/result.png" alt="Drawing" style="width: 500px;"/>

Will produce *lot/result.png*
Requires *plot/solution.dat* to work properly.

## Other
Use:
- `make clean` to clean up
- `make frames [dim=%d] [itr=%d] [frames=%d]` to create *images/\*.png* frames using serial version
- `make gif [delay=%d]` to create GIF, requires *images/\*.png* frames to work properly

Where `frames` number of images produced and `delay` time delay between frames in the GIF.

<img src="./aux/animation.gif" alt="Drawing" style="width: 500px;"/>
