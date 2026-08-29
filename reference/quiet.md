# Suppress stdout, messages, and warnings

Evaluates an expression while suppressing both stdout and message
output. This is necessary because some dependencies (like RAS) emit
verbose output directly to stderr or through R's message system.

## Usage

``` r
quiet(expr)
```

## Arguments

- expr:

  Expression to evaluate.

## Value

The result of evaluating `expr`.
