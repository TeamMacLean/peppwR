# Create a new peppwr_power object

Create a new peppwr_power object

## Usage

``` r
new_peppwr_power(mode, question, answer, simulations, params, call)
```

## Arguments

- mode:

  Either "aggregate" or "per_peptide"

- question:

  What was solved for: "power", "sample_size", or "effect_size"

- answer:

  The computed answer

- simulations:

  List of simulation details

- params:

  List of input parameters used

- call:

  Original function call

## Value

A peppwr_power object
