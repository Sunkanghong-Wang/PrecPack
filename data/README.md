# Benchmark instances

This directory contains the shared item instances and labeled graphs distributed with PrecPack. It contains no PrecPack result rows, solution files, experiment records, logs, or timing data.

## Mathematical interpretation

Each item or task $i$ has a positive size $w_i$ and is assigned to a zero-based bin or station position $b_i$. Every position has capacity $C$. For each precedence arc $(i,j,t_{ij})$, feasibility requires

$$
b_j-b_i \ge t_{ij}.
$$

PrecPack minimizes the number of positions from position 0 through the last used position. The three supported problems differ only in the arc labels:

| Problem | Separation labels |
| --- | --- |
| SALBP-I | $t_{ij}=0$: predecessor and successor may share a station when their task order is feasible. |
| BPP-P | $t_{ij}=1$: every successor must be placed in a strictly later bin. |
| BPP-GP | Nonnegative labels are read from a `.graph` file. |

## Directory map and counts

| Path | Contents | Count |
| --- | --- | ---: |
| `instances/otto/n_XXXX/` | Shared Otto base instances for $n\in\{20,50,100,250,500,750,1000\}$ | 3,675 |
| `instances/otto/n_0050_permuted/` | Nine additional task-order permutations of each Otto-50 instance | 4,725 |
| `instances/scholl/` | Shared classical Scholl instances | 269 |
| `bpp-gp-graphs/separation-01/n_XXXX/` | Labeled graphs with separations in $\{0,1\}$ | 3,675 |
| `bpp-gp-graphs/separation-03/n_XXXX/` | Labeled graphs with separations in $\{0,1,2,3\}$ | 3,675 |

The `.txt` files are shared across problem types. BPP-P uses the seven Otto instance-size sets and the Scholl set. SALBP-I uses those same files together with `n_0050_permuted`. BPP-GP pairs the Otto base files with one of the two labeled-graph sets under `bpp-gp-graphs/`. Problem semantics are selected by `--problem`; no symbolic links or duplicate instance files are required.

## Instance `.txt` format

Instance files use the established ALB text format and the lowercase `.txt` suffix. They are UTF-8 or ASCII text. Blank lines and lines beginning with `#` are ignored. Task identifiers are one-based.

```text
<number of tasks>
4

<cycle time>
10

<order strength>
0.5

<task times>
1 4
2 3
3 5
4 2

<precedence relations>
1,3
2,4

<end>
```

Required sections are:

- `<number of tasks>`: one positive integer $n$;
- `<cycle time>`: one positive capacity $C$;
- `<task times>`: exactly one `task_id weight` row for every task, with positive weights;
- `<precedence relations>`: zero or more `from,to` rows defining direct arcs.

`<order strength>` is optional metadata. In SALBP-I every ALB arc receives label 0. In BPP-P every ALB arc receives label 1. When a BPP-GP graph is supplied, its arc set replaces the ALB precedence section.

## Labeled `.graph` format

Graph files are ASCII text. The first nonblank line is the exact number of arcs. Every remaining nonblank row contains three integers:

```text
3
1 3 0
2 4 1
3 4 2
```

The columns are `from`, `to`, and nonnegative separation $t_{ij}$. Item identifiers are one-based. The declared count must equal the number of arc rows. Negative labels, missing labels, invalid item identifiers, and cyclic graphs are rejected.

## Pairing BPP-GP files

An instance file and a labeled graph pair when they have the same size directory and stem. For example:

```text
instances/otto/n_0020/instance_n=20_1.txt
bpp-gp-graphs/separation-03/n_0020/instance_n=20_1.graph
```

## Sources and acknowledgements

- The Scholl instances come from the public [1993 SALBP benchmark set](https://assembly-line-balancing.de/salbp/benchmark-data-sets-1993/).
- The Otto SALBP instances come from the public [2013 SALBP benchmark set](https://assembly-line-balancing.de/salbp/benchmark-data-sets-2013/) accompanying Otto, Otto, and Scholl, [“Systematic data generation and test design for solution algorithms on the example of SALBPGen for assembly line balancing”](https://doi.org/10.1016/j.ejor.2012.12.029), *European Journal of Operational Research* 228(1), 33–45, 2013.
- The BPP-P and BPP-GP benchmark sets follow the construction described by Kramer, Dell'Amico, and Iori, [“A batching-move iterated local search algorithm for the bin packing problem with generalized precedence constraints”](https://doi.org/10.1080/00207543.2017.1341065), *International Journal of Production Research* 55(21), 6288–6304, 2017.

We gratefully thank the authors and maintainers of these benchmark sets for making them available to the research community and for supporting reproducible research in packing and assembly-line optimization.
