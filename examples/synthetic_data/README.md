# Overview

This folder contains very simple examples, mainly for testing but also as explanations.

One example is the column AAAACCCC, together with trees of different kinds, for example the balanced tree
(((a,b),(c,d)),((e,f),(g,h))), and compared with the completely unbalanced tree
(a,(b,(c,(d,(e,(f,(g,h))))))). What logo do you expect?

Also try ACCCCCC with the two trees. Here, PIC should make the unbalanced tree produce a completely different visualization than the balanced one, right?


# Data

## ex1

*   The first two columns are completely conserved: only "A".
*   Then come two columns with a full spread of "ARVESTID".
*   Columns 5 and 6 contain half A and half C in the upper and lower halves, respectively.
*   Columns 7 and 8 alternate between A and C on every other row.
*   Columns 9 and 10 have "R" everywhere except in the first (a) and last (h) sequence.
    -   In this case, the outcome should differ greatly depending on which tree is used.


### Trees

* `ex1_t1.tree`: Fully balanced tree
* `ex1_t2.tree`: Completely unbalanced, with "a" as its own subtree under the root.
* `ex1_t3.tree`: Completely unbalanced, like `ex1_t2`, but with "h" instead of "a" 
   as sole subtree under the root.


```
$ for x in 1 2 3; do nw_display ex1_t$x".tree"; done
                                                    +------------------------+ a
                          +-------------------------+
                          |                         +------------------------+ b
 +------------------------+
 |                        |                         +------------------------+ c
 |                        +-------------------------+
 |                                                  +------------------------+ d
=|
 |                                                  +------------------------+ e
 |                        +-------------------------+
 |                        |                         +------------------------+ f
 +------------------------+
                          |                         +------------------------+ g
                          +-------------------------+
                                                    +------------------------+ h

 |------------|-----------|------------|------------|-----------|-------------
 0         0.05         0.1         0.15          0.2        0.25
 substitutions/site

 +---------------------------------------------------------------------------+ a
 |
=|          +----------------------------------------------------------------+ b
 |          |
 +----------+          +-----------------------------------------------------+ c
            |          |
            +----------+          +------------------------------------------+ d
                       |          |
                       +----------+         +--------------------------------+ e
                                  |         |
                                  +---------+          +---------------------+ f
                                            |          |
                                            +----------+          +----------+ g
                                                       +----------+
                                                                  +----------+ h

 |---------------------|--------------------|---------------------|-----------
 0                   0.1                  0.2                   0.3
 substitutions/site

 +---------------------------------------------------------------------------+ h
 |
=|          +----------------------------------------------------------------+ b
 |          |
 +----------+          +-----------------------------------------------------+ c
            |          |
            +----------+          +------------------------------------------+ d
                       |          |
                       +----------+         +--------------------------------+ e
                                  |         |
                                  +---------+          +---------------------+ f
                                            |          |
                                            +----------+          +----------+ g
                                                       +----------+
                                                                  +----------+ a

 |---------------------|--------------------|---------------------|-----------
 0                   0.1                  0.2                   0.3
```
