# readHAC

R package to read acoustic HAC format

## Installation

To install `readHAC` run
```
devtools::install_github("kaskr/HAC",subdir="readHAC")
```

## Tested tuples

| Tuple | Tested | TODO               |
|-------|--------|--------------------|
|    20 |   TRUE |      NA            |
|    41 |   TRUE |      NA            |
|    42 |  FALSE |      NA            |
|   100 |  FALSE |      NA            |
|   200 |  FALSE |      NA            |
|   210 |   TRUE |      NA            |
|   901 |   TRUE |      NA            |
|  1000 |  FALSE |      NA            |
|  1001 |  FALSE |      NA            |
|  2000 |  FALSE |      NA            |
|  2001 |  FALSE |      NA            |
|  2002 |  FALSE |      NA            |
|  2100 |   TRUE |      NA            |
|  4000 |   TRUE |      NA            |
|  9001 |   TRUE |      NA            |
| 10000 |   TRUE |      NA            |
| 10001 |   TRUE |      NA            |
| 10010 |  FALSE |      NA            |
| 10011 |  FALSE | split RLE samples  |
| 10030 |   TRUE |      NA            |
| 10031 |  FALSE |      NA            |
| 10040 |   TRUE |      NA            |
| 10090 |   TRUE |      NA            |
| 10100 |  FALSE |      NA            |
| 10140 |   TRUE |      NA            |
| 10142 |  FALSE |      NA            |
| 11000 |  FALSE |      NA            |
| 65397 |  FALSE |      NA            |
| 65534 |   TRUE |      NA            |
| 65535 |   TRUE |      NA            |
